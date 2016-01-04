from collections import defaultdict
import peaks
import peak_groups
import chrom
import numpy as np
import savitzky_golay as sg
import parameters

__author__ = 'Tiannan Guo, ETH Zurich 2015'


class NestedDict(defaultdict):
    def __init__(self):
        super(NestedDict, self).__init__(NestedDict)

    def __deepcopy__(self):
        return self


class Chromatogram(object):

    def __init__(self, rt_list_three_values_csv, i_list_csv):

        rt_list = map(float, peaks.rt_three_values_to_full_list_string(rt_list_three_values_csv).split(','))
        i_list = map(float, i_list_csv.split(','))

        # i_list_smoothed = smooth_chromatogram_using_Savitzky_Golay(i_list)

        self.rt_list = rt_list
        self.i_list = i_list
        # self.i_list_smoothed = i_list_smoothed

        max_peaks, __ = peaks.peakdetect(i_list, rt_list, 2.0, 0.3)
        # max_peaks_smoothed, __ = peaks.peakdetect(i_list_smoothed, rt_list, 6.0, 0.3)

        if len(max_peaks) > 0:

            # max_peaks_smoothed = filter_smoothed_peaks_based_on_raw_peaks(max_peaks, max_peaks_smoothed)

            max_peaks_all = filter_peaks_based_on_peak_shape(max_peaks, i_list, rt_list)

            self.peak_apex_rt_list = [rt for (rt, i) in max_peaks_all]
            self.peak_apex_i_list = [i for (rt, i) in max_peaks_all]

        else:
            # if no peak found. Most likely there is no signal. Use looser criteria to detect peaks

            max_peaks, __ = peaks.peakdetect(i_list, rt_list, 1.0, 0.3)
            # max_peaks_smoothed, __ = peaks.peakdetect(i_list_smoothed, rt_list, 1.0, 0.3)

            # max_peaks_smoothed = filter_smoothed_peaks_based_on_raw_peaks(max_peaks, max_peaks_smoothed)

            max_peaks_all = filter_peaks_based_on_peak_shape(max_peaks, i_list, rt_list)

            if len(max_peaks) > 0:
                self.peak_apex_rt_list = [rt for (rt, i) in max_peaks_all]
                self.peak_apex_i_list = [i for (rt, i) in max_peaks_all]

            else:
                # if still no peak found, most likely it is an empty chrom.
                # use the point in the median value as peak
                self.peak_apex_rt_list = [np.median(np.array(rt_list))]
                self.peak_apex_i_list = [np.median(np.array(i_list))]


def filter_peaks_based_on_peak_shape(max_peaks, i_list, rt_list):

    max_peaks_all = []

    max_peaks_all = filter_peaks_based_on_peak_shape_worker(max_peaks, i_list, rt_list, max_peaks_all)

    if len(max_peaks_all) == 0:
        # if no peak is found, find the best peak in the chrom
        max_peaks_all = filter_peaks_based_on_peak_shape_worker2(max_peaks, i_list, rt_list, max_peaks_all)

    return max_peaks_all


def filter_peaks_based_on_peak_shape_worker2(max_peaks, i_list, rt_list, max_peaks_all):

    max_fold_change_score = -1
    best_rt = -1
    best_i = -1

    for rt, i in max_peaks:

        rt_left, rt_right = chrom.get_peak_boundary(rt_list, i_list, rt)
        i_apex = float(i)
        i_left = peak_groups.get_intensity_for_closest_rt(rt_left, rt_list, i_list)
        i_right = peak_groups.get_intensity_for_closest_rt(rt_right, rt_list, i_list)
        fold_change_left = i_apex / (i_left + 1.0)
        fold_change_right = i_apex / (i_right + 1.0)

        fold_change_score = fold_change_left * fold_change_right

        if fold_change_score > max_fold_change_score:
            max_fold_change_score = fold_change_score
            best_rt = rt
            best_i = i

    max_peaks_all.append((best_rt, best_i))

    return max_peaks_all


def filter_peaks_based_on_peak_shape_worker(max_peaks, i_list, rt_list, max_peaks_all):

    for rt, i in max_peaks:

        rt_left, rt_right = chrom.get_peak_boundary(rt_list, i_list, rt)
        i_apex = float(i)
        i_left = peak_groups.get_intensity_for_closest_rt(rt_left, rt_list, i_list)
        i_right = peak_groups.get_intensity_for_closest_rt(rt_right, rt_list, i_list)
        fold_change_left = i_apex / (i_left + 1.0)
        fold_change_right = i_apex / (i_right + 1.0)

        if fold_change_left >= parameters.PEAK_SHAPE_FOLD_VARIATION_CRUDE and fold_change_right >= parameters.PEAK_SHAPE_FOLD_VARIATION_CRUDE:
            max_peaks_all.append((rt, i))

    return max_peaks_all


def filter_smoothed_peaks_based_on_raw_peaks(peaks, peaks_smoothed):

    peaks_smoothed2 = []

    for rt, i in peaks_smoothed:
        if_select = decide_whether_choose_a_smoothed_rt(rt, [rt for (rt, i) in peaks])
        if if_select == 1:
            peaks_smoothed2.append((rt, i))

    return peaks_smoothed2


def decide_whether_choose_a_smoothed_rt(rt0, rt_list):

    if_select = 0

    rt_closest_left = find_closest_rt_left(rt0, rt_list)
    rt_closest_right = find_closest_rt_right(rt0, rt_list)

    if abs(rt_closest_left - rt0) > 10 and abs(rt_closest_right - rt0) > 10:
        if_select = 1
    return if_select


def find_closest_rt_left(rt0, rt_list):

    rt1 = rt_list[0]

    rt_dif = 999.999

    for rt in rt_list[1:]:
        if rt < rt0:
            rt_dif2 = rt0 - rt
            if rt_dif2 < rt_dif:
                rt1 = rt
                rt_dif = rt_dif2
    return rt1


def find_closest_rt_right(rt0, rt_list):
    rt1 = rt_list[0]
    rt_dif = 999.999
    for rt in rt_list[1:]:
        if rt > rt0:
            rt_dif2 = rt - rt0
            if rt_dif2 < rt_dif:
                rt1 = rt
                rt_dif = rt_dif2
    return rt1


def smooth_chromatogram_using_Savitzky_Golay(i_list):

    i_list2 = sg.savitzky_golay(np.array(i_list), 7, 3)  # window size 11, polynomial order 3. Optimized for chrom
    return i_list2.tolist()


class ReferenceSample(object):

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


class PeakGroup(object):

    def __init__(self, chrom_data, tg, sample, rt):
        self.rt = rt
        (self.matched_fragments, self.matched_fragments_rt,
         self.matched_fragments_i,
         self.matched_fragments_peak_rt_left,
         self.matched_fragments_peak_rt_right) = peak_groups.find_matched_fragments(chrom_data, tg, sample, rt)

        self.num_matched_fragments = len(self.matched_fragments)
        self.if_ms1_peak = peak_groups.check_if_ms1_peak(chrom_data, tg, sample, rt)
