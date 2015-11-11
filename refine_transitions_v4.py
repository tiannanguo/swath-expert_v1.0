__author__ = 'Tiannan Guo, ETH Zurich 2015'

import csv
import peakdetect
import numpy as np
import gzip
import sys
import time
from collections import defaultdict
class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

PEAK_TOLERANCE = 6  #seconds, use to find if the peak is found

id_mapping_file = 'goldenSets90.txt'
in_file = 'allChrom_1.txt.gz'
out_file =in_file.repalce('allChrom', 'goodChrom')

# in_file = sys.argv [1] #"test_chromGz_com120.txt.gz"
# out_file = "refined_" + in_file #  "test_chromGz_com120_refined.txt.gz"


def read_id_file():
    id = []
    id_mapping = {}
    with open(id_mapping_file) as i:
        reader = csv.reader(i, delimiter="\t")
        for row in reader:
            id.append(row[1])
            id_mapping[row[0].lower()] = row[1]
    return id, id_mapping

def refine_transitions(id, id_mapping):
    d = ngram()
    with gzip.open(in_file, 'rb') as i:
        reader = csv.DictReader(i, delimiter="\t")
        for row in reader:


def main():
    # if the in_file is empty, do nothing
    num_lines = gzip.open(in_file).read().count('\n')
    if (num_lines <= 1):
        pass
    else:
        id, id_mapping = read_id_file()
        d = refine_transitions(id, id_mapping)


start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))


# transition_groups = []
transition_list ={}
sample_name = {}

with gzip.open(in_file, 'rb') as IN_FILE, gzip.open (out_file, 'wb') as OUT_FILE:
    reader = csv.DictReader (IN_FILE, delimiter ="\t")
    writer = csv.DictWriter(OUT_FILE, fieldnames=reader.fieldnames, delimiter = "\t")
    headers = {}
    for n in writer.fieldnames:
        headers[n] = n
    writer.writerow(headers)

    sample_name = [s[:-3] for s in reader.fieldnames if "_SW_rt" in s]
    # print "samples include ", sample_name

    for row in reader:

        if (row['transition_name'] == row['transition_group_id']):
            writer.writerow(row)
        else:
            best_sample = row['best_sample']
            best_rt= float (row['best_rt'])

            rt_column_name = best_sample + "_rt"
            int_column_name = best_sample + "_int"

            rt_list = map(float, row[rt_column_name].split(","))

            rt_list2 = list (np.arange(rt_list[0], rt_list[1],rt_list[2]))
            if (rt_list[1] - rt_list2[-1] >1):
                rt_list2.append(rt_list[1])
            int_list = map(float, row[int_column_name].split(","))

            peak_max = []
            peak_min = []
            if (len(rt_list2) != len(int_list)):
                print "ERROR in transition name is %s " % transition_name
                print "rt_list length %d is different from int_list length %d" % (len(rt_list2), len(int_list))
                print "rt_list is %s" % rt_list
                print "rt_list2 is %s" % rt_list2
                print "int_list is %s" % int_list

            else:
                peak_max, peak_min = peakdetect.peakdetect(int_list, rt_list2, 10.0, 0.3)

            for peak_max_pair in peak_max:
                # print "type of peak_max_rt is ", type(peak_max_pair[0])
                if abs(peak_max_pair[0] - best_rt) < PEAK_TOLERANCE:
                    #a peak is found for the transition in the rt range, write it to the outfile
                    writer.writerow(row)
                    break

