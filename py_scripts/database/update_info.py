#!/usr/bin/python
from peewee import *
import sys
sys.path.insert(0, "/home/anna/bioinformatics/ngs/py_scripts/")
from blast.classes.blast_parser import BlastParser
from database.models import *
from common_helpers.parse_csv import *

def update_function_info(csv_info_path):
    info_dict = csv_to_dict(csv_info_path, main_key='seqid')
    with db.atomic():
        for seqid in info_dict:
            seq_info = info_dict[seqid]
            other_new_info = {}
            for key in seq_info:
                if key not in ('function', 'mitoscore'):
                    other_new_info[key] = seq_info[key]
            sequence = Sequence.get(Sequence.seqid == seqid)
            extra_data = sequence.extra_data
            extra_data.update(other_new_info)
            query = Sequence.update(function=seq_info['function'], mitoscore=seq_info['mitoscore'], extra_data=extra_data).where(Sequence.seqid==seqid)
            query.execute()
    return 0

def update_len_ratio():
    i = 0
    for hit in BlastHit.select():
        extra_data = hit.extra_data
        length = hit.length
        hit.alen_qlen = float(length/float(extra_data['qlen']))
        hit.alen_slen = float(length/float(extra_data['slen']))
        hit.save()
    return 0


csv_info_path = '/home/anna/bioinformatics/euglenozoa/mitocarta/Human.MitoCarta.2.0.csv'

update_len_ratio()