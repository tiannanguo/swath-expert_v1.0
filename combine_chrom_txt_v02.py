__author__ = 'Tiannan Guo, ETH Zurich 2015'

import gzip
import csv
import time
import sys
from collections import defaultdict
class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

id_mapping_file = sys.argv[1]
window_num = sys.argv[2]

# id_mapping_file = "goldenSets90_test.txt"
# window_num = 10

title = ['transition_name','transition_group_id', 'best_rt', 'best_sample', 'best_score', 'protein', 'irt', 'precursor_mz', 'transition_mz', 'transition_i']

def read_id_file():
    id = []
    id_mapping = {}
    with open(id_mapping_file) as i:
        reader = csv.reader(i, delimiter="\t")
        for row in reader:
            id.append(row[1])
            id_mapping[row[0].lower()] = row[1]
    return id, id_mapping

def read_in_ms1_ms2_data(i, d, id_mapping):
    in_file = i + '_' + str(window_num) + '.chrom.txt.gz'
    tg = []
    with gzip.open(in_file, 'rb') as gz:
        reader = csv.DictReader(gz, delimiter="\t")
        for row in reader:
            transition = row[title[0]]
            tg.append(row[title[1]])
            for j in range(1,10):
                d[transition][title[j]] =row[title[j]]
            if (d[transition]['best_sample'].lower() in id_mapping):
                d[transition]['best_sample'] = id_mapping[d[transition]['best_sample'].lower()]
            d[transition]['rt_list'][i] = row['rt_list']
            d[transition]['i_list'][i] = row['i_list']

    in_file2 = i + '_ms1.chrom.txt.gz'
    with gzip.open(in_file2, 'rb') as gz:
        reader = csv.DictReader(gz, delimiter="\t")
        for row in reader:
            transition = row[title[1]]
            if (transition in tg):
                for j in range(1,10):
                    d[transition][title[j]] = row[title[j]]
                if (d[transition]['best_sample'].lower() in id_mapping):
                    d[transition]['best_sample'] = id_mapping[d[transition]['best_sample'].lower()]
                d[transition]['rt_list'][i] = row['rt_list']
                d[transition]['i_list'][i] = row['i_list']

    return d

def print_data(id, d, window_num):
    out_file = 'com_chrom_' + str(window_num) + '.txt.gz'
    with gzip.open(out_file, 'wb') as o:
        writer = csv.writer(o, delimiter="\t")
        for e in id:
            title.append(e + '_rt')
            title.append(e + '_i')

        writer.writerow(title)

        for transition in d.keys():
            dd = []
            dd.append(transition)
            for j in range(1,10):
                dd.append(d[transition][title[j]])
            for e in id:
                if e in d[transition]['rt_list']:
                    dd.append(d[transition]['rt_list'][e])
                    dd.append(d[transition]['i_list'][e])
                else:
                    dd.append('NA')
                    dd.append('NA')

            writer.writerow(dd)

def main():
    id, id_mapping = read_id_file()
    d = ngram()

    for i in id:
        d = read_in_ms1_ms2_data(i, d, id_mapping)

    print_data(id, d, window_num)

start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))