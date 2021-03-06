#!python3
from Bio import Entrez
Entrez.email = "a.nenarokova@gmail.com"
from encoder import XML2Dict

def read_diamond_result(diamond_taxa_result_path):
    taxid_dict = {}
    with open(diamond_taxa_result_path) as diamond_f:
        for line in diamond_f:
            line_split = line.rstrip().split("\t")
            qseqid = line_split[outfmt_opt_list.index("qseqid")]
            taxid = line_split[outfmt_opt_list.index("taxid")]
            evalue = float(line_split[outfmt_opt_list.index("evalue")])
            if taxid in taxid_dict:
                taxid_dict[taxid]["count"] += 1
            else:
                taxid_dict[taxid] = {}
                taxid_dict[taxid]["count"] = 1
    return taxid_dict

diamond_taxa_result_path = "/home/anna/bioinformatics/myxo_local/myxo/smol_diamond_nr_tax2.tsv"
outfmt_opts = "qseqid taxid evalue"
outfmt_opt_list = outfmt_opts.split(" ")

outpath = "/home/anna/bioinformatics/myxo_local/myxo/tax_analysis_result.tsv"

taxid_dict = read_diamond_result(diamond_taxa_result_path)

taxids = list(taxid_dict.keys())
# print (taxids)

ncbi_xml_record = Entrez.efetch(db="taxonomy", id=taxids).read()

ncbi_dict = XML2Dict().parse(ncbi_xml_record)

taxa_dict_list = ncbi_dict['TaxaSet']['Taxon']

for taxon in taxa_dict_list:
    taxid = taxon['TaxId'].decode("utf-8")
    name = taxon['ScientificName'].decode("utf-8")
    try:
        lineage = taxon['Lineage'].decode("utf-8")
    except:
        lineage = taxon['Lineage']
    taxid_dict[taxid]["name"] = name
    taxid_dict[taxid]["lineage"] = lineage

taxid_dict['0']["name"] = "none"
taxid_dict['0']["lineage"] = "none"

with open(outpath, "w") as outfile:
    header = "taxid\tcount\tname\tlineage\n"
    outfile.write(header)
    for taxid in taxid_dict:
        count = taxid_dict[taxid]["count"]
        name = taxid_dict[taxid]["name"]
        lineage = taxid_dict[taxid]["lineage"]
        row = "{}\t{}\t{}\t{}\n".format(taxid,count,name,lineage)
        outfile.write(row)






