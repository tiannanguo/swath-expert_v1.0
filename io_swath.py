__author__ = 'guot'

import numpy as np
from collections import defaultdict
import csv
import gzip
import chrom
class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

MIN_TRANSITION_NUMBER = 5  # after optimization, I find set min 5 transitions is good


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
    d = ngram()
    with gzip.open(chrom_file, 'rb') as i:
        r = csv.DictReader(i, delimiter='\t')
        for row in r:
            fragment = row[title[0]]
            tg = row[title[1]]
            d[tg]['reference_sample']['peak_rt'] = float(row['best_rt'])
            d[tg]['reference_sample']['name'] = row['best_sample']
            d[tg]['reference_sample']['score'] = float(row['best_score'])
            d[tg][title[5]] = row[title[5]]
            d[tg][title[6]] = float(row[title[6]])
            d[tg][title[7]] = float(row[title[7]])
            if fragment == tg: #ms1
                for k in sample_id:
                    d[tg]['precursor']['rt_list'][k] = row[k + '_rt']
                    d[tg]['precursor']['i_list'][k] = row[k + '_i']
            else: #ms2
                for j in range(8, 10):
                    d[tg]['fragments'][fragment][title[j]] = float(row[title[j]])
                for k in sample_id:
                    d[tg]['fragments'][fragment]['rt_list'][k] = row[k + '_rt']
                    d[tg]['fragments'][fragment]['i_list'][k] = row[k + '_i']
    return d








