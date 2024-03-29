#!/usr/bin/python
import shutil
import os

def listdir_nohidden(path):
	for f in os.listdir(path):
		if not f.startswith('.'):
			yield f

def read_list(list_path):
	result_list = []
	with open (list_path) as list_file:
		for line in list_file:
			result_list.append(line.rstrip())
	return result_list

def copy_files(ids_list, indir, outdir, delimiter="."):
	for file_name in listdir_nohidden(indir):
		id = file_name.split(delimiter)[0]
		if id in ids_list:
			print ("copying", file_name)
			file_path = indir + file_name
			outpath = outdir + file_name
			shutil.copyfile(file_path, outpath)
	return 0

keep_list = ["EP00002","EP00004","EP00023","EP00026","EP00029","EP00031","EP00033","EP00039","EP00057","EP00067","EP00072","EP00083","EP00098","EP00103","EP00104","EP00107","EP00110","EP00119","EP00120","EP00127","EP00129","EP00131","EP00134","EP00145","EP00150","EP00157","EP00159","EP00162","EP00163","EP00164","EP00167","EP00176","EP00185","EP00202","EP00206","EP00210","EP00229","EP00247","EP00248","EP00260","EP00269","EP00271","EP00279","EP00298","EP00300","EP00314","EP00333","EP00379","EP00390","EP00397","EP00454","EP00455","EP00457","EP00466","EP00473","EP00503","EP00513","EP00609","EP00611","EP00656","EP00667","EP00671","EP00681","EP00686","EP00696","EP00697","EP00698","EP00727","EP00736","EP00738","EP00741","EP00742","EP00749","EP00750","EP00753","EP00759","EP00761","EP00762","EP00770","EP00771","EP00794","EP00797","EP00802","EP00808","EP00810","EP00813","EP00820","EP00852","EP00857","EP00865","EP00895","EP00900","EP00924","EP00927","EP00931","EP01027","EP01028","EP01062","EP01080","EP01082","EP01083","EP01084","EP01086","EP01087","EP01091","EP01094","EP01104","EP01113","EP01114","EP01117","EP01125","EP01127","EP01130","EP01134","EP01146"]
indir = "/Users/vl18625/work/euk/protein_sets/eukprot/eukprot3/proteins/"
outdir = "/Users/vl18625/work/euk/protein_sets/anna_dataset/anna_eukprot3_v3_21_06_23/"
copy_files(keep_list, indir, outdir, delimiter="_")
