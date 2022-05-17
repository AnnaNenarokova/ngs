#!/usr/bin/python3
import re
import os
from Bio import SeqIO

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

def parse_hmmreport(hmm_report_path, columns_str=False):
	results = []
	if not columns_str:
		columns_str = "sseqid s_accession qseqid q_accession evalue score bias 1_domain_evalue 1_domain_score 1_domain_bias exp reg clu ov env dom rep inc"
		columns_list = columns_str.split(" ")
	with open(hmm_report_path) as hmmfile:
		for line in hmmfile:
			if line[0] != "#":
				line_edited = re.sub(' +', ' ', line)
				line_split = line_edited.rstrip().split(" ")
				qseqid = line_split[columns_list.index("qseqid")]
				sseqid = line_split[columns_list.index("sseqid")]
				evalue = float(line_split[columns_list.index("evalue")])
				bitscore = float(line_split[columns_list.index("score")])
				result_dict = {"qseqid": qseqid, "sseqid": sseqid, "evalue": evalue, "bitscore": bitscore}
				results.append(result_dict)
	return results

def prepare_hmm_dict(hmm_report_dir, hmm_ext, proteome_ext, n_best, max_evalue, monobranch=False):
	hmm_dict = {}
	for hmm_report in listdir_nohidden(hmm_report_dir):
		hmm_report_name_split = hmm_report.split(proteome_ext)
		proteome_file = hmm_report_name_split[0] + proteome_ext
		hmm_file = hmm_report_name_split[1].split(hmm_ext)[0]
		if monobranch:
			hmm_file_split = hmm_file.split("_")
			cog_file = hmm_file_split[0] + ("_") + hmm_file_split[1] + ".faa"
		else:
			cog_file = hmm_file
		if proteome_file not in hmm_dict:
			hmm_dict[proteome_file] = {}

		if cog_file not in hmm_dict[proteome_file].keys():
			hmm_dict[proteome_file][cog_file] = {}
		hmm_report_path = hmm_report_dir + hmm_report
		hmm_results = parse_hmmreport(hmm_report_path)
		hmm_results = list(filter(lambda result: (result["evalue"] <= max_evalue), hmm_results))
		if len(hmm_results) > n_best:
			hmm_results = sorted(hmm_results, key=lambda result: result["evalue"])
			hmm_results = hmm_results[:n_best]
		for hmm_result in hmm_results:
			hmm_dict[proteome_file][cog_file][hmm_result["sseqid"]] = ""
	return hmm_dict

def prepare_hmm_dict_old(hmm_report_dir, hmm_ext, proteome_ext, n_best, max_evalue):
	hmm_dict = {}
	for hmm_report in listdir_nohidden(hmm_report_dir):
		hmm_report_name_split = hmm_report.split(proteome_ext)
		proteome_file = hmm_report_name_split[0] + proteome_ext
		cog_file = hmm_report_name_split[1].split(hmm_ext)[0]
		if proteome_file not in hmm_dict:
			hmm_dict[proteome_file] = {}
		hmm_dict[proteome_file][cog_file] = {}
		hmm_report_path = hmm_report_dir + hmm_report
		hmm_results = parse_hmmreport(hmm_report_path)
		hmm_results = list(filter(lambda result: (result["evalue"] <= max_evalue), hmm_results))
		if len(hmm_results) > n_best:
			hmm_results = sorted(hmm_results, key=lambda result: result["evalue"])
			hmm_results = hmm_results[:n_best]
		for hmm_result in hmm_results:
			hmm_dict[proteome_file][cog_file][hmm_result["sseqid"]] = ""
	return hmm_dict

def prepare_fastas(hmm_report_dir, proteome_dir, cog_dir, result_dir, hmm_ext=".hmm", proteome_ext=".fasta", n_best=1, max_evalue=0.00001):
	print("Parcing hmm_reports")
	hmm_dict = prepare_hmm_dict(hmm_report_dir, n_best=n_best, max_evalue=max_evalue)
	print("Parcing proteomes sequences")
	for proteome in listdir_nohidden(proteome_dir):
		proteome_dict = hmm_dict[proteome]
		proteome_path = proteome_dir + proteome
		for record in SeqIO.parse(proteome_path, "fasta"):
			sseqid = record.id
			for cog in proteome_dict:
				if record.id in proteome_dict[cog]:
					hmm_dict[proteome][cog][sseqid] = record
	print("Writing down sequences")
	for cog_file in listdir_nohidden(cog_dir):
		cog = cog_file.split(".faa")[0]
		out_records = []
		cog_path = cog_dir + cog_file
		outpath = result_dir + cog_file
		for record in SeqIO.parse(cog_path, "fasta"):
			out_records.append(record)
		for proteome in hmm_dict:
			proteome_cog_dict = hmm_dict[proteome][cog]
			for sseqid in proteome_cog_dict:
				out_records.append(proteome_cog_dict[sseqid])
		SeqIO.write(out_records, outpath, "fasta")
	return 0

def prepare_fastas_keep_list(hmm_report_dir, proteome_dir, cog_dir, result_dir, monobranch=False, keep_list_path=False, hmm_ext=".txt", proteome_ext=".fasta", n_best=1, max_evalue=0.00001):
	if keep_list_path: keep_list = read_list(keep_list_path)
	print("Parcing hmm_reports")
	hmm_dict = prepare_hmm_dict(hmm_report_dir, hmm_ext, proteome_ext, n_best, max_evalue, monobranch=monobranch)
	print("Parcing proteome sequences")
	proteome_set = set()
	for proteome in listdir_nohidden(proteome_dir):
		proteome_id = proteome.split("-")[0]
		if not keep_list_path or (proteome_id in keep_list):
			proteome_set.add(proteome)
			proteome_dict = hmm_dict[proteome]
			proteome_path = proteome_dir + proteome
			for record in SeqIO.parse(proteome_path, "fasta"):
				sseqid = record.id
				for cog_file in proteome_dict:
					if record.id in proteome_dict[cog_file]:
						hmm_dict[proteome][cog_file][sseqid] = record
	print("Writing down sequences")
	for cog_file in listdir_nohidden(cog_dir):
		out_set = set()
		out_records = []
		cog_path = cog_dir + cog_file
		outpath = result_dir + cog_file
		for record in SeqIO.parse(cog_path, "fasta"):
			seqid = record.id
			if seqid not in out_set:
				out_records.append(record)
				out_set.update(seqid)
		for proteome in proteome_set:
			if cog_file in hmm_dict[proteome]:
				proteome_cog_dict = hmm_dict[proteome][cog_file]
				for sseqid in proteome_cog_dict:
					record = proteome_cog_dict[sseqid]
					seqid = record.id
					if seqid not in out_set:
						out_records.append(record)
						out_set.update(seqid)
		SeqIO.write(out_records, outpath, "fasta")
	return 0


hmm_report_dir = "/scratch/nenarokova/euk/markers/be_mono_results/alpha/hmm_results/"
proteome_dir = "/scratch/nenarokova/euk/proteomes/anna_set_16_05_22/"
cog_dir = "/scratch/nenarokova/euk/markers/bacteria/BacEuk_markers_faa/"
result_dir = "/scratch/nenarokova/euk/markers/be_mono_results/alpha/faa/"

prepare_fastas_keep_list(hmm_report_dir, proteome_dir, cog_dir, result_dir)
