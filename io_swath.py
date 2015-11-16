__author__ = 'guot'

import numpy as np
from collections import defaultdict
import csv
import gzip
import chrom
import parameters
import peaks

class Nested_dict(defaultdict):
    def __init__(self):
        super(Nested_dict, self).__init__(Nested_dict)

class Chromatogram(object):

    def __init__(self, rt_list_three_values_csv, i_list_csv):

        rt_list = map(float, peaks.rt_three_values_to_full_list_string(rt_list_three_values_csv).split(','))
        i_list = map(float, i_list_csv.split(','))

        self.rt_list = rt_list
        self.i_list = i_list

        max_peaks, __ = peaks.peakdetect(i_list, rt_list, 9.0, 0.3)

        self.peak_apex_rt = [rt for (rt, i) in max_peaks]
        self.peak_apex_i = [i for (rt, i) in max_peaks]

    def size(self):
        return len(self.rt_list)


class Reference_sample(object):

    def __init__(self, sample_name, score, peak_rt):
        self.sample_name = sample_name
        self.score = score
        self.peak_rt = peak_rt
        self.peak_rt_left = ''
        self.peak_rt_right = ''

    def get_peak_boundary(self, peak_rt_left, peak_rt_right):
        self.peak_rt_left = peak_rt_left
        self.peak_rt_right = peak_rt_right

class Peptide(object):

    def __init__(self, protein_name, irt, ms1_dict, ms2_dict):
        self.protein = protein_name
        self.irt = irt
        self.ms1_dict = {}
        self.ms2_dict = {}

    def get_ms1_data(self, ms1_name, ms1_mz):
        self.ms1_dict[ms1_name] = ms1_mz

    def get_ms2_data(self, ms2_name, ms2_mz):
        self.ms1_dict[ms2_name] = ms2_mz

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
    chrom_data = Nested_dict()
    peptide_data = Nested_dict()
    with gzip.open(chrom_file, 'rb') as i:
        r = csv.DictReader(i, delimiter='\t')
        for row in r:
            fragment = row['transition_name']
            tg = row['transition_group_id']
            peak_rt = float(row['best_rt'])
            ref_sample_name = row['best_sample']
            ref_sample_score = float(row['best_score'])

            ref_sample_data[tg] = Reference_sample(ref_sample_name, ref_sample_score, peak_rt)

            peptide_data[tg]['ms1']['name'] = tg
            peptide_data[tg]['ms1']['mz'] = float(row['precursor_mz'])
            peptide_data[tg]['ms1']['protein'] = row['protein']
            peptide_data[tg]['ms1']['protein'] = float(row['irt'])

            peptide_data[tg]['ms2']['name'].append(fragment)
            peptide_data[tg]['ms2']['mz'].append(float(row['transition_mz']))

            for k in sample_id:
                rt_list_three_values_csv = row[k + '_rt']
                i_list_csv = row[k + '_i']
                chrom_data[tg][k][fragment] = Chromatogram(rt_list_three_values_csv, i_list_csv)

    return ref_sample_data, chrom_data, peptide_data








