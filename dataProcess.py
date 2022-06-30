import os
import sqlite3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def encodeNuc(nuc):
    enc = ["a", "c", "g", "t", ".", "-"]
    i = enc.index(nuc)
    one_hot = [0] * 6
    one_hot[i] = 1
    return one_hot


def getSNPdf(db_name):
    # get tSNE input data from snp database
    # establish connection to sqlite database
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # get all possible sequence_names aka types
    sql_command = "SELECT DISTINCT sequence_name FROM unambiguous ORDER BY sequence_name;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    sequence_names = list()
    for seq in content:
        sequence_names.append(seq[0])

    # get all reference_GenWidePos for SNPs
    sql_command = "SELECT DISTINCT reference_GenWidePos FROM unambiguous ORDER BY reference_GenWidePos;"
    cursor.execute(sql_command)
    content = cursor.fetchall()
    positions = list()
    for pos in content:
        positions.append(pos[0])

    query = []
    for sequence in sequence_names:
        # query SNP pattern and position from database
        sql_command = "SELECT sequence_name, SNP_pattern, reference_GenWidePos FROM unambiguous WHERE" \
                      " unambiguous.sequence_name='" + str(sequence) + "' ORDER BY reference_GenWidePos;"
        cursor.execute(sql_command)
        content = cursor.fetchall()
        query.append(content)
    flat_query = [x for xs in query for x in xs]
    df = pd.DataFrame(flat_query, columns=['Sequence_name', 'SNP', 'Position'])
    df['Reference'] = df['SNP'].str[1]
    df['SNP'] = df['SNP'].str[0]
    return df


def encDF(df, default_enc, common_pos):
    symb_count = len(default_enc)
    seq_count = df['Sequence_name'].nunique()
    sequences = df['Sequence_name'].unique().tolist()

    # create training vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    reference_3D = np.full((1, len(common_pos), symb_count), default_enc)
    vec_3D = np.full((seq_count, len(common_pos), symb_count), default_enc)

    for index, row in df.iterrows():
        if row['Position'] in common_pos:
            seq = sequences.index(row["Sequence_name"])
            pos = common_pos.index(row["Position"])

            enc_entry = encodeNuc(row["SNP"])
            vec_3D[seq][pos] = enc_entry

            reference_3D[0][pos] = encodeNuc(row["Reference"])

    vec_3D = np.append(vec_3D, reference_3D, axis=0)

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = vec_3D.shape
    vec_2D = vec_3D.reshape((sample, position * one_hot_enc))
    return vec_2D


def encQuery(df, default_enc, common_pos):
    sample_names = list(df.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    # create vector query SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    eval_3D = np.full((len(sample_names), len(common_pos), len(default_enc)), default_enc)

    s = 0
    for sample in sample_names:
        p = 0
        for entry in df[sample]:
            pos = common_pos.index(df['Position'].iloc[p])
            eval_3D[s][pos] = encodeNuc(entry.lower())
            p += 1
        s += 1

    # reshape training vector from 3D to 2D
    sample, position, one_hot_enc = eval_3D.shape
    eval_2D = eval_3D.reshape((sample, position * one_hot_enc))
    return eval_2D


def assignStrains(sample_names):
    # all samples manually grouped by strain
    group1 = ['BosniaA']
    group2 = ['Fribourg', 'CDC2', 'GHA1', 'Gauthier', 'IND1', 'SAM1', 'SamoaD']
    group3 = ['Seattle81', 'SEA86', 'NE20', 'BAL3', 'BAL73', 'Chicago', 'NIC1', 'NIC2', 'Dallas']

    strain_labels = ['TEN', 'TPE', 'Nicols (TPA)', 'SS14 (TPA)']

    strains = []
    for sample in sample_names:
        if sample in group1:
            strain = strain_labels[0]
        elif sample in group2:
            strain = strain_labels[1]
        elif sample in group3:
            strain = strain_labels[2]
        else:
            strain = strain_labels[3]
        strains.append(strain)
    return strains


def getADAdataset(tsvFile, db_name, default_enc):
    # get all SNPs from reference database
    SNPdf = getSNPdf(db_name)
    ref_positions = SNPdf["Position"].unique().tolist()
    ref_positions.sort()

    # read tsvFile containing query sequence SNPs
    if os.path.isfile(tsvFile):
        query = pd.read_csv(tsvFile, sep='\t', header=1)
    # create vector of SNPs, with field for each SNP site (reference_GenWidePos)
    sample_names = list(query.columns)
    sample_names.remove('Position')
    sample_names.remove('Reference')

    # find common SNP positions in reference database and query
    query_pos = query.Position.tolist()
    common_pos = list(set(map(str, ref_positions)).intersection(query_pos))
    common_pos.sort(key=lambda pos: int(pos))

    # remove all rows from dataframes that don't have common reference genome position for SNP
    query = query[query['Position'].isin(common_pos)]
    query['Position'] = pd.to_numeric(query['Position'])
    query.reset_index()
    common_pos = list(map(int, common_pos))
    SNPdf = SNPdf[SNPdf['Position'].isin(common_pos)]
    SNPdf.reset_index()

    sequence_names = SNPdf['Sequence_name'].unique().tolist()
    sequence_names.sort()
    sequence_names.append("Reference(NC_021490)")

    # create training and test vector of SNPs, with field for each SNP site (reference_GenWidePos) using one-hot encoding
    train_2D = encDF(SNPdf, default_enc, common_pos)
    test_2D = encQuery(query, default_enc, common_pos)
    dataset = {'train': train_2D, 'train_seq_names': sequence_names, 'test': test_2D, 'test_sample_names': sample_names}
    return dataset

# com_pos = list(map(int, ['6442', '7179', '8831', '9319', '9560', '9984', '18609', '20497', '20773', '41132', '47527', '59894', '64080', '64081', '71972', '73158', '79676', '79965', '85612', '88008', '88684', '94075', '94901', '101869', '122630', '132121', '132281', '134933', '134957', '134965', '134974', '134975', '134981', '134990', '134994', '134996', '135001', '135002', '135008', '135246', '135250', '135975', '135996', '135997', '136003', '136007', '139557', '148115', '149166', '149360', '149537', '150557', '151017', '151207', '151218', '151220', '152311', '152511', '152620', '152648', '152652', '152661', '153158', '153200', '153206', '153222', '153227', '153233', '153236', '155812', '156721', '157504', '158149', '158163', '158167', '158177', '158178', '158180', '158189', '158193', '158201', '158202', '158271', '158274', '158275', '158339', '158340', '158341', '158346', '158348', '158354', '158363', '158364', '158384', '158386', '158415', '158426', '158590', '158618', '158619', '158624', '158631', '158649', '158714', '158721', '158774', '158787', '158796', '158882', '158884', '158885', '158915', '158916', '158976', '159084', '159147', '159156', '159251', '159297', '159308', '159309', '159312', '159320', '159323', '159374', '159376', '159499', '159650', '169148', '170012', '171515', '174137', '187064', '187073', '187177', '197405', '198040', '198145', '198186', '198194', '198269', '198334', '198352', '198396', '198428', '202251', '202410', '213742', '226275', '235204', '237746', '237860', '240069', '242666', '246989', '258741', '264881', '267834', '271107', '271108', '278740', '279291', '283649', '287971', '288872', '295415', '295416', '302121', '302720', '312085', '312087', '314485', '319090', '319233', '319251', '319354', '319648', '320692', '320834', '320912', '320968', '327063', '327426', '327602', '327666', '330866', '331040', '331044', '331045', '331048', '331504', '331512', '332760', '332781', '332782', '332788', '332792', '333343', '333351', '333357', '333551', '333876', '333877', '333879', '333936', '333939', '333940', '333944', '334118', '334483', '336221', '336910', '338109', '342664', '347027', '347031', '347034', '347036', '347061', '347062', '347068', '347075', '347083', '347260', '347401', '347807', '347822', '347824', '347825', '347827', '347928', '347955', '347956', '353352', '353649', '354963', '357127', '364849', '370050', '370681', '370688', '372549', '372899', '373221', '374505', '378996', '388147', '408751', '413262', '414832', '415693', '417189', '420034', '421764', '422205', '424072', '428595', '435853', '442656', '449518', '459976', '460450', '462617', '462742', '462743', '462982', '462983', '463004', '463007', '463011', '463042', '463043', '463425', '471391', '477634', '492583', '492772', '492773', '492826', '492835', '492836', '492842', '492845', '492850', '492869', '492895', '492905', '492924', '492974', '492975', '492979', '492982', '492994', '492995', '492999', '493004', '493013', '493030', '493057', '493060', '493061', '493112', '493204', '493216', '493237', '493285', '493288', '493308', '493358', '493374', '493394', '493531', '493549', '493577', '493582', '493601', '493602', '493605', '493614', '495695', '500171', '501035', '503421', '512349', '514619', '514720', '514823', '514885', '514894', '514897', '515326', '517755', '521090', '521150', '521754', '522981', '523014', '523015', '523017', '523031', '523042', '523060', '523146', '523203', '523291', '523308', '523330', '523380', '523395', '523468', '523471', '523480', '523485', '523488', '523506', '523531', '523545', '523612', '523620', '523632', '523669', '523723', '524113', '528814', '530139', '537571', '537846', '538516', '539783', '544427', '550331', '555872', '556275', '556289', '556803', '556820', '556971', '557589', '557598', '557603', '557612', '557618', '557619', '557670', '557768', '557771', '564544', '573412', '574564', '576927', '577258', '579791', '586227', '592290', '592614', '593153', '593424', '593429', '593435', '593436', '593438', '593439', '593441', '593444', '593452', '593477', '593563', '593686', '593738', '593740', '593742', '593747', '593750', '593753', '593754', '594031', '594032', '594033', '594035', '594038', '594040', '594044', '594138', '594143', '594144', '594184', '594187', '594215', '594221', '594222', '594571', '594930', '602216', '606165', '606171', '606229', '606552', '606591', '606939', '607169', '608437', '613421', '617495', '619959', '622191', '623981', '624844', '626141', '629984', '630004', '630013', '630020', '636378', '642894', '648150', '652210', '656469', '656856', '657057', '661257', '661263', '661299', '661366', '663195', '663307', '664641', '667340', '671127', '672262', '672264', '672269', '672324', '673752', '673773', '673774', '673780', '673784', '674335', '674343', '674349', '674542', '675109', '675143', '675148', '675152', '675272', '675274', '675548', '675559', '675560', '675653', '675656', '685611', '690421', '690671', '703109', '703137', '708563', '710654', '729217', '730935', '737426', '737724', '741934', '749757', '760212', '767985', '772966', '773215', '773571', '773572', '783503', '789700', '790531', '790985', '799933', '800591', '803821', '804726', '809311', '811875', '812238', '816454', '821958', '825071', '831568', '837713', '858819', '860410', '861564', '870687', '870734', '874149', '887590', '890154', '906101', '910995', '912782', '916216', '928188', '929715', '935052', '935053', '935096', '935102', '935554', '936048', '936177', '936178', '936230', '936235', '936241', '936242', '936250', '936256', '936259', '936431', '936440', '936443', '936450', '936473', '936529', '936613', '936637', '936642', '936648', '936654', '936655', '936658', '936725', '936771', '936774', '936781', '936788', '936794', '936823', '936829', '936841', '936851', '936862', '936880', '936892', '937067', '937149', '937164', '937180', '945224', '945292', '945308', '945331', '945520', '945526', '945527', '945541', '945542', '945693', '945696', '945702', '945703', '945705', '945707', '945709', '945830', '945838', '945845', '945851', '945853', '946103', '946112', '946118', '946126', '946127', '946138', '946142', '946298', '948473', '948666', '951392', '951595', '955916', '956280', '960746', '964252', '970792', '975837', '975840', '975883', '975920', '976105', '976106', '976114', '976115', '976116', '976118', '976119', '976120', '976121', '976122', '976123', '976124', '976125', '976126', '976127', '976128', '976131', '976243', '976244', '976246', '976247', '976248', '976249', '976250', '976254', '976255', '976256', '976257', '976258', '976259', '976260', '976261', '976263', '976264', '976265', '976266', '976267', '976269', '976271', '976275', '976379', '976380', '976381', '976382', '976383', '976384', '976390', '976391', '976392', '976393', '976396', '976404', '976407', '976408', '976410', '976411', '976412', '976413', '976416', '976422', '976425', '976510', '976512', '976514', '976522', '976529', '976585', '976588', '976589', '976590', '976591', '976592', '976595', '976597', '976598', '976601', '976603', '976604', '976612', '976616', '976620', '976623', '976629', '976647', '976684', '976692', '976726', '976728', '976729', '976732', '976733', '976735', '976741', '976744', '976746', '976747', '976748', '976753', '976754', '976755', '976756', '976759', '976760', '977220', '977221', '977812', '977900', '978458', '979718', '985097', '992054', '992764', '992765', '994358', '997714', '1006195', '1006197', '1007869', '1023862', '1027228', '1035113', '1042556', '1049545', '1049618', '1049623', '1050211', '1050223', '1050442', '1050539', '1051257', '1051898', '1051901', '1051902', '1051944', '1052027', '1052069', '1052230', '1052231', '1052302', '1052414', '1052671', '1052755', '1052919', '1052923', '1052929', '1052974', '1053298', '1053366', '1053547', '1053617', '1054009', '1057409', '1057766', '1059461', '1059494', '1059629', '1062042', '1063680', '1072433', '1079844', '1081597', '1086484', '1102208', '1103911', '1113171', '1113176', '1124233', '1125380', '1125514', '1125555', '1125741', '1125865', '1125866', '1125870', '1125877', '1125924', '1126002', '1126028', '1126720', '1126729', '1126730', '1126748', '1126753', '1126771', '1126783', '1126784', '1126786', '1126789', '1126793', '1126795', '1126813', '1126819', '1126823', '1126828', '1126846', '1126847', '1127003', '1127023', '1127024', '1127162', '1127185', '1127186', '1127189', '1127190', '1127201', '1127202']))
# print(encDF(getSNPdf("snps.db"), [0, 0, 0, 0, 1, 0], com_pos))
