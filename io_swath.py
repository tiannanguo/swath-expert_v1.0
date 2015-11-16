__author__ = 'guot'

import numpy as np
from collections import defaultdict
import csv
import gzip
import chrom
import parameters
import peaks
import data_holder



def read_id_file(id_mapping_file):
    sample_id = []
    id_mapping = {}
    with open(id_mapping_file) as i:
        reader = csv.reader(i, delimiter="\t")
        for row in reader:
            sample_id.append(row[1])
            id_mapping[row[0].lower()] = row[1]
    return sample_id, id_mapping

def read_com_chrom_file(chrom_file, sample_id, title):
    ref_sample_data = {}
    chrom_data = data_holder.Nested_dict()
    peptide_data = data_holder.Nested_dict()
    with gzip.open(chrom_file, 'rb') as i:
        r = csv.DictReader(i, delimiter='\t')
        for row in r:
            fragment = row['transition_name']
            tg = row['transition_group_id']
            peak_rt = float(row['best_rt'])
            ref_sample_name = row['best_sample']
            ref_sample_score = float(row['best_score'])

            ref_sample_data[tg] = data_holder.Reference_sample(ref_sample_name, ref_sample_score, peak_rt)

            peptide_data[tg]['ms1']['preMz'] = float(row['precursor_mz'])
            peptide_data[tg]['ms1']['protein'] = row['protein']
            peptide_data[tg]['ms1']['irt'] = float(row['irt'])

            peptide_data[tg]['ms2'][fragment] = float(row['transition_mz'])

            for k in sample_id:
                rt_list_three_values_csv = row[k + '_rt']
                i_list_csv = row[k + '_i']
                chrom_data[tg][k][fragment] = data_holder.Chromatogram(rt_list_three_values_csv, i_list_csv)

    return ref_sample_data, chrom_data, peptide_data








