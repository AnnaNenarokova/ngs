#!/usr/bin/python3
from bioservices.kegg import KEGG

output = open('/home/kika/MEGAsync/blasto_project/genes/replication/lmaj.txt', 'w')

kegg = KEGG()
# pathway = kegg.get('tbr03030')
pathway = kegg.get('lma03030')
dict_data = kegg.parse(pathway)
# print(dict_data)

genes = {}
for key in dict_data.keys():
	if key == 'GENE':
		genes = dict_data[key]

for key in genes.keys():
	output.write(key + '\n')



# for key in dict_data['ORTHOLOGY'].keys():
# 	value = dict_data['ORTHOLOGY'][key]
# 	if 'EC' in value:
# 		EC = value.split('[EC:')[1]
# 		EC = EC.split(']')[0]
# 		EC = EC.replace(' ', '\n')
# 		output.write(EC + '\n')
# output.close()

# # g = x.get('tbr03440:Tb11.01.0910/aaseq')
# # print(g)

# # res = x.parse_kgml_pathway("tbr03440")
# # print(res['entries'][0])

# # for key, value in dict_data['GENE'].items():
# # 	print(key, value)

# # for gene in dict_data['GENE']:
# # 	output.write(gene + '\n')

# kegg = KEGG()
# pathway = kegg.get('ath00564')
# dict_data = kegg.parse(pathway)

# for value in dict_data['GENE'].values():
# 	EC = value.split('[EC:')[1]
# 	EC = EC.replace(']', '')
# 	EC = EC.replace(' ', '\n')
# 	output.write(EC + '\n')