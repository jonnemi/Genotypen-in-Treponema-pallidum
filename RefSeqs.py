from Bio import Entrez
Entrez.email = 'jonathan.nemitz@student.uni-tuebingen.de'

handle = Entrez.esearch(db="nucleotide", retmax=100,term="Treponema pallidum pallidum AND srcdb_refseq[PROP]", idtype="acc")
record = Entrez.read(handle)
idList = record["IdList"]
toShort = ["NR_076531.1", "NR_076156.1", "NR_075513.1", "NR_075106.1"]
for id in toShort:
    idList.remove(id)


savePath = r"C:/Users/Jonathan/Documents/BioInformatik/6. Semester/Bachelorarbeit/Mauve_alignments/mauveIn/"

def fetch_id(path, id):
    filename = path + id + ".gbk"
    with Entrez.efetch(db="nuccore", rettype="gb", style="withparts", retmode="text", id=id) as handle:
        out_handle = open(filename, "w")
        out_handle.write(handle.read())
        out_handle.close()
        handle.close()

def fetch_all():
    reference = "NC_021490.2"
    if reference in record["IdList"]:
        refPath = "C:/Users/Jonathan/Documents/BioInformatik/6. Semester/Bachelorarbeit/Mauve_alignments/mauveIn/Reference_"
        fetch_id(refPath, reference)
        record["IdList"].remove(reference)

    count = 1
    for id in idList:
        fetch_id(savePath, id)
        total = len(idList)
        print("Fetch " + str(count) + " of " + str(total) + " done.")
        count += 1

print(record)