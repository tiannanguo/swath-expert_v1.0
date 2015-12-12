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

def read_com_chrom_file(chrom_file, sample_id):
    ref_sample_data = {}
    chrom_data = data_holder.Nested_dict()
    peptide_data = data_holder.Nested_dict()

    # sometimes, in the golden standard data set, multiple "best_sample"s are found. The best in water may not be the best in human
    # in this case, check the input file and write to a new file
    # select the best sample with lowest m_score
    chrom_file2 = use_one_best_sample(chrom_file)

    with gzip.open(chrom_file2, 'rb') as i:
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

                ###########
                if k == 'gold80' and fragment.startswith('1191_'):
                    pass
                chrom_data[tg][k][fragment] = data_holder.Chromatogram(rt_list_three_values_csv, i_list_csv)

    return ref_sample_data, chrom_data, peptide_data


def get_tg_list(chrom_file):
    tg_list = {}
    with gzip.open(chrom_file, 'rb') as i:
        r = csv.DictReader(i, delimiter="\t")
        for row in r:
            tg_list[row['transition_group_id']] = 1
    return tg_list.keys()

def get_best_sample_for_each_tg(chrom_file, tg_list):

    best_sample = {}
    best_score = {}
    best_rt = {}

    for tg in tg_list:
        best_sample[tg] = ''
        best_score[tg] = 1.0
        best_rt[tg] = -1.0

        with gzip.open(chrom_file, 'rb') as i:

            r = csv.DictReader(i, delimiter="\t")

            for row in r:
                tg2 = row['transition_group_id']
                if tg2 == tg and float(row['best_score']) < float(best_score[tg]):
                    best_sample[tg] = row['best_sample']
                    best_rt[tg] = row['best_rt']
                    best_score[tg] = float(row['best_score'])

    return best_sample, best_score, best_rt

def use_one_best_sample(chrom_file):

    chrom_file2 = chrom_file.replace('.txt.gz', '_2.txt.gz')

    # get all tg
    tg_list = get_tg_list(chrom_file)

    # get the best sample for each tg
    best_sample, best_score, best_rt = get_best_sample_for_each_tg(chrom_file, tg_list)

    with gzip.open(chrom_file, 'rb') as i, gzip.open(chrom_file2, 'wb') as o:
        r = csv.DictReader(i, delimiter="\t")
        fieldnames1 = r.fieldnames
        w = csv.DictWriter(o, delimiter="\t", fieldnames=fieldnames1)
        w.writeheader()
        # write data
        for row in r:
            tg = row['transition_group_id']
            row['best_sample'] = best_sample[tg]
            row['best_rt'] = best_rt[tg]
            row['best_score'] = best_score[tg]
            w.writerow(row)

    return chrom_file2







