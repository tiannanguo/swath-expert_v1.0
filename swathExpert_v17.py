__author__ = 'Tiannan Guo, ETH Zurich 2015'

import gzip
import sys
import csv
import numpy as np
import time
from collections import defaultdict
class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

import io_swath
import peaks
import peak_groups
import r_code
import chrom

MIN_FRAGMENTS = 4
MAX_RT_TOLERANCE = 10
MAX_PEAK_WIDTH = 100
PEAK_WIDTH_FOLD_VARIATION = 2.0
PEAK_SHAPE_FOLD_VARIATION = 3.0

# name of files
chrom_file = 'allChrom_1.txt.gz'    #sys.argv[1]
id_mapping_file = 'goldenSets90.txt'
out_file = chrom_file.repalce('.txt.gz', '.R')
out_file_poor_tg = chrom_file.replace('.txt.gz', '.poor.txt')
quant_file = chrom_file.replace('.txt.gz', '.quant.txt')

# parameters
PEAK_TOLERANCE = 6  # seconds, use to find if the peak is found
title = ['transition_name','transition_group_id', 'best_rt', 'best_sample', 'best_score', 'protein', 'irt', 'precursor_mz', 'transition_mz', 'transition_i']



def main():

    # read input file of sample inforamtion
    sample_id, id_mapping = io_swath.read_id_file()

    # read input chrom file, construct the d class
    # write out rt_list, i_list, etc
    d = io_swath.read_com_chrom_file(chrom_file, sample_id, title)

    # based on rt_list and i_list, find peaks for all chrom, write into d class: "fragments", "precursor"
    d = peaks.find_peaks(d, sample_id)

    # based on peaks of fragments, find out all peak groups, write into d class: "peak_groups"
    d = peak_groups.find_peak_groups(d, sample_id, MIN_FRAGMENTS, MAX_RT_TOLERANCE)

    # based on peak groups found in the reference sample, find out fragments that form good peaks, remove the rest fragments
    d = peak_groups.refine_fragments_based_on_reference_sample(d, MAX_RT_TOLERANCE, MIN_FRAGMENTS)

    # compute the peak boundary for the reference sample, write to display_pg
    d = chrom.compute_reference_sample_peak_boundary(d, MAX_PEAK_WIDTH)

    # based on the display_pg for reference sample, get the best matched peak groups from all other samples
    # and then write to display_pg
    d = peak_groups.find_best_peak_group_based_on_reference_sample(d, sample_id, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)

    # compute peak area for display_pg
    d = chrom.compute_peak_area_for_all(d)

    # write r code
    all_r_codes = r_code.write_r_code_for_all_samples(d, sample_id)


start_time = time.time()
main()
print "--- %s seconds ---" % (time.time() - start_time)
