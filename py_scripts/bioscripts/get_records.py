#!/usr/bin/python
from Bio import SeqIO

l = [
"EG_transcript_8890",
"EG_transcript_9775",
"EG_transcript_20525",
"EG_transcript_20372",
"EG_transcript_13787",
"EG_transcript_10977",
"EG_transcript_8499",
"EG_transcript_14580",
"EG_transcript_10058",
"EG_transcript_26432",
"EG_transcript_10781",
"EG_transcript_14903",
"EG_transcript_18984",
"EG_transcript_24306",
"EG_transcript_33335",
"EG_transcript_43502",
"EG_transcript_12354",
"EG_transcript_13951",
"EG_transcript_14568",
"EG_transcript_13173",
"EG_transcript_17919",
"EG_transcript_37915",
"EG_transcript_2635",
"EG_transcript_14436",
"EG_transcript_9541",
"EG_transcript_14883",
"EG_transcript_3952",
"EG_transcript_12688",
"EG_transcript_16636",
"EG_transcript_21190",
"EG_transcript_20549",
"EG_transcript_32878",
"EG_transcript_27754",
"EG_transcript_13739",
"EG_transcript_18872",
"EG_transcript_20704",
"EG_transcript_43870",
"EG_transcript_32586",
"EG_transcript_13361",
"EG_transcript_24179",
"EG_transcript_70244",
"EG_transcript_20494",
"EG_transcript_12447",
"EG_transcript_8645",
"EG_transcript_13449",
"EG_transcript_14243",
"EG_transcript_6376",
"EG_transcript_11471",
"EG_transcript_14864",
"EG_transcript_11639",
"EG_transcript_18956",
"EG_transcript_11556",
"EG_transcript_13179",
"EG_transcript_12961",
"EG_transcript_34430",
"EG_transcript_12294",
"EG_transcript_46692",
"EG_transcript_66034",
"EG_transcript_43713",
"EG_transcript_49602",
"EG_transcript_36144",
"EG_transcript_57265",
"EG_transcript_60368",
"EG_transcript_38927",
"EG_transcript_11705",
"EG_transcript_12866",
"EG_transcript_33292",
"EG_transcript_16863",
"EG_transcript_14142",
"EG_transcript_16247",
"EG_transcript_7118",
"EG_transcript_14662",
"EG_transcript_10950",
"EG_transcript_17253",
"EG_transcript_13456",
"EG_transcript_36501",
"EG_transcript_24421",
"EG_transcript_7568",
"EG_transcript_39807",
"EG_transcript_14563",
"EG_transcript_13498",
"EG_transcript_11072",
"EG_transcript_49841",
"EG_transcript_12754",
"EG_transcript_8382",
"EG_transcript_14285",
"EG_transcript_28264",
"EG_transcript_58517",
"EG_transcript_52164",
"EG_transcript_23496",
"EG_transcript_62041",
"EG_transcript_39609",
"EG_transcript_38110",
"EG_transcript_32981",
"EG_transcript_22203",
"EG_transcript_22943",
"EG_transcript_15359",
"EG_transcript_18724",
"EG_transcript_19661",
"EG_transcript_15192",
"EG_transcript_22946",
"EG_transcript_13542",
"EG_transcript_9059",
"EG_transcript_6915",
"EG_transcript_17405",
"EG_transcript_14562",
"EG_transcript_11449",
"EG_transcript_12949",
"EG_transcript_9296",
"EG_transcript_23394",
"EG_transcript_39893",
"EG_transcript_17905",
"EG_transcript_9818",
"EG_transcript_22436",
"EG_transcript_26453",
"EG_transcript_26301",
"EG_transcript_19816",
"EG_transcript_23268",
"EG_transcript_21588",
"EG_transcript_9443",
"EG_transcript_24487",
"EG_transcript_47186",
"EG_transcript_12142",
"EG_transcript_14008",
"EG_transcript_8533",
"EG_transcript_11591",
"EG_transcript_13467",
"EG_transcript_38272",
"EG_transcript_8248",
"EG_transcript_23999",
"EG_transcript_29447",
"EG_transcript_23123",
"EG_transcript_13873",
"EG_transcript_17331",
"EG_transcript_18758",
"EG_transcript_33770",
"EG_transcript_72005",
"EG_transcript_28440",
"EG_transcript_17931",
"EG_transcript_62987",
"EG_transcript_69335",
"EG_transcript_22374",
"EG_transcript_59404",
"EG_transcript_39706",
"EG_transcript_28515",
"EG_transcript_13593",
"EG_transcript_26735",
"EG_transcript_16435",
"EG_transcript_49389",
"EG_transcript_53971",
"EG_transcript_46569",
"EG_transcript_27927",
"EG_transcript_45419",
"EG_transcript_41848",
"EG_transcript_29771",
"EG_transcript_28505",
"EG_transcript_41201",
"EG_transcript_41841",
"EG_transcript_19975",
"EG_transcript_69249",
"EG_transcript_67124",
"EG_transcript_57639",
"EG_transcript_65067",
"EG_transcript_58459",
"EG_transcript_22648",
"EG_transcript_33101",
"EG_transcript_41744",
"EG_transcript_69930",
"EG_transcript_69571",
"EG_transcript_46805",
"EG_transcript_46809",
"EG_transcript_32320",
"EG_transcript_48855",
"EG_transcript_2555",
"EG_transcript_51012",
"EG_transcript_25297",
"EG_transcript_18745",
"EG_transcript_29816",
"EG_transcript_32467",
"EG_transcript_19369",
"EG_transcript_20521",
"EG_transcript_24577",
"EG_transcript_16112",
"EG_transcript_38034",
"EG_transcript_34278",
"EG_transcript_26056",
"EG_transcript_31675",
"EG_transcript_20436",
"EG_transcript_50571",
"EG_transcript_28065",
"EG_transcript_21020",
"EG_transcript_27241"

]

fasta = '/home/anna/Dropbox/PhD/bioinformatics/euglena/data/euglena_all_proteins.fasta'
results = []

for record in SeqIO.parse(fasta, "fasta"):
    if record.id in l:
        results.append(record)

outpath = '/home/anna/Dropbox/PhD/bioinformatics/euglena/data/MTERF_euglena.faa'

SeqIO.write(results, outpath, "fasta")
