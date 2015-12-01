from collections import defaultdict
import peaks
import peak_groups
import chrom
import numpy as np
import copy

__author__ = 'Tiannan Guo, ETH Zurich 2015'

class Nested_dict(defaultdict):
    def __init__(self):
        super(Nested_dict, self).__init__(Nested_dict)

    def __deepcopy__(self):
        return self

class Chromatogram(object):

    def __init__(self, rt_list_three_values_csv, i_list_csv):

        rt_list = map(float, peaks.rt_three_values_to_full_list_string(rt_list_three_values_csv).split(','))
        i_list = map(float, i_list_csv.split(','))

        self.rt_list = rt_list
        self.i_list = i_list

        max_peaks, __ = peaks.peakdetect(i_list, rt_list, 8.0, 0.3)

        if len(max_peaks) > 0:

            self.peak_apex_rt_list = [rt for (rt, i) in max_peaks]
            self.peak_apex_i_list = [i for (rt, i) in max_peaks]

        else:
            # if no peak found. Most likely there is no signal. Use looser criteria to detect peaks

            max_peaks, __ = peaks.peakdetect(i_list, rt_list, 1.0, 0.3)

            if len(max_peaks) > 0:
                self.peak_apex_rt_list = [rt for (rt, i) in max_peaks]
                self.peak_apex_i_list = [i for (rt, i) in max_peaks]

            else:
                # if still no peak found, most likely it is an empty chrom.
                # use the point in the median value as peak
                self.peak_apex_rt_list = [np.median(np.array(rt_list))]
                self.peak_apex_i_list = [np.median(np.array(i_list))]

        # peak_rt_left, peak_rt_right = peaks.get_all_peak_boundary(rt_list, i_list, max_peaks)

        # self.peak_apex_rt = peak_rt_left
        # self.peak_apex_right = peak_rt_right



    # def size(self):
    #     return len(self.rt_list)


class Reference_sample(object):

    def __init__(self, sample_name, score, peak_rt):
        self.sample_name = sample_name
        self.score = score
        self.peak_rt = peak_rt
        self.peak_rt_found = ''
        self.peak_rt_left = ''
        self.peak_rt_right = ''

    def read_peak_boundary(self, peak_rt_left, peak_rt_right):
        self.peak_rt_left = peak_rt_left
        self.peak_rt_right = peak_rt_right

    def read_peak_rt_found(self, rt_found):
        self.peak_rt_found = rt_found


class Peak_group(object):
    def __init__(self, chrom_data, tg, sample, rt, fold_change):
        self.rt = rt
        self.matched_fragments, self.matched_fragments_rt, self.matched_fragments_i,\
            self.matched_fragments_peak_rt_left, self.matched_fragments_peak_rt_right = \
                peak_groups.find_matched_fragments(chrom_data, tg, sample, rt, fold_change)

        self.num_matched_fragments = len(self.matched_fragments)
        self.if_ms1_peak = peak_groups.check_if_ms1_peak(chrom_data, tg, sample, rt)


