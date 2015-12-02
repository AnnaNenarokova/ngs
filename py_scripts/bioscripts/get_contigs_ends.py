from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
f = ''
start_len = 130 
for record in SeqIO.parse(f, "fasta"):
	seq = record.seq
	record_len = len(seq)
	description = 'length ' + str(record_len)
	if record_len <= 2*end_len:
		seq_out = SeqRecord(seq, id = str(record.id), description = description)
		out.append(seq_out)
	else:
		seq1_out = SeqRecord(seq[0:end_len], id =  str(record.id) + ' ' + 'first ' + str(end_len), description = description)
		out.append(seq1_out)
		seq2_out = SeqRecord(seq[-end_len:-1], id = str(record.id) + ' ' + 'last ' + str(end_len),  description = description)
		out.append(seq2_out)
SeqIO.write(out, file_out, "fasta")