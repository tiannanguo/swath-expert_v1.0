__author__ = 'Tiannan Guo, ETH Zurich 2015'


import numpy as np
from collections import defaultdict
import csv
import gzip
import math
import peaks


def swath_peak_group_filter_expert(in_file, transition_groups, this_tg, sample_id,
                                   reference_sample_peak_group_info):

    # get the list of all peak groups
    rt_list, int_list = file_io.read_rt_int_list_for_other_samples(in_file, this_tg, sample_id)

    peaks_rt = peakdetect.find_peaks_for_a_sample(rt_list, int_list)

    peak_groups = process_peak_groups.find_peak_groups(peaks_rt)

    # expert system to filter out poor peak groups
    num_filtered_peak_group, peak_groups_filtered = \
        swath_peak_group_filter_expert_worker(this_tg,
                                              rt_list,
                                              int_list,
                                              peak_groups,
                                              len(transition_groups[this_tg]),
                                              reference_sample_peak_group_info)

    peak_info = {}

    # one good peak group is found
    if num_filtered_peak_group == 1 and peak_groups_filtered.keys()[0] > 0:
        # NOTE: peak_groups_filtered.keys()[0] > 0, because if there is only no peak, the defaultDict still takes 0.0 as the first value,
        # in this case, there is no peak group

        # compute peak boundary
        peak_info = compute_peak_boundary_and_area_expert(rt_list, int_list, transition_groups, this_tg,
                                                          peak_groups_filtered.keys()[0], sample_id)

    # if no peak group left, use the peak from the reference sample. RULE x
    elif num_filtered_peak_group <= 1:
        peak_info['sample_id'] = sample_id
        peak_info['peak_group_rt_left'] = reference_sample_peak_group_info['peak_group_rt_left']
        peak_info['peak_group_rt_apex'] = reference_sample_peak_group_info['peak_group_rt_apex']
        peak_info['peak_group_rt_right'] = reference_sample_peak_group_info['peak_group_rt_right']
        peak_info['peak_group_area'] = 0.0

    else:
        print "ERROR: swath_peak_group_filter_expert: sth is wrong. More than 1 peak groups passed the filter"
################
      # some bug????

    return num_filtered_peak_group, peak_info


def swath_peak_group_filter_expert_worker(this_tg, rt_list, int_list, peak_group,
                                          reference_sample_num_transitions,
                                          reference_sample_peak_group_info):

    # select the best peak group, afterwards, perform a quality check

    num_good_transitions, peak_group2 = count_number_of_good_transition_expert(
        peak_group, reference_sample_num_transitions)

    peak_group_filtered = defaultdict(list)

    # RULE 1: if no peak group does not have enough transitions, return the empty dict
    if num_good_transitions == 0:
        peak_group_filtered = peak_group2

    else:

        number_highly_correlated_peak_groups, peak_group3 = compute_fragment_ranking_expert(this_tg, rt_list, int_list,
                                                                                            peak_group2, reference_sample_peak_group_info)

        # RULE 2: check fragment ranking
        if number_highly_correlated_peak_groups == 0:
            peak_group_filtered = peak_group3

        elif number_highly_correlated_peak_groups == 1:  # one good peak group
            peak_group_filtered = peak_group3

        else:
            # RULE 3: check MS1
            num_good_ms1, peak_group4 = examine_ms1_expert(this_tg, peak_group3)

            if num_good_ms1 == 1:
                peak_group_filtered = peak_group4

            # if there is no ms1 or multiple ms1, then examine peak shape, select the
            # one with highest score
            elif num_good_ms1 == 0 or num_good_ms1 > 1:

                # RULE 4: if multiple peaks found, eg. 76_AAALEAMK_2 , sample MCF7_a,
                # select the one with better peak shape
                peak_group5 = find_best_peak_shape_expert(
                    this_tg, rt_list, int_list, peak_group4, reference_sample_peak_group_info)

                peak_group_filtered = peak_group5

    num_filtered_peak_group = len(peak_group_filtered.keys())

    return num_filtered_peak_group, peak_group_filtered


def count_number_of_good_transition_expert(peak_group, best_sample_num_transitions):

    num_good_transition = 0

    peak_group2 = defaultdict(list)

    # get the peak groups with highest numbers of transitions
    max_transition = 0
    for rt in peak_group:
        this_peak_group_transition_number = len(peak_group[rt]) - 1
        if this_peak_group_transition_number > max_transition:
            max_transition = this_peak_group_transition_number

    min_transition = max_transition - 3

    if min_transition < 5:  # must at least 5 transitions
        min_transition = 5

    for rt in peak_group:
        this_peak_group_transition_number = len(peak_group[rt]) - 1
        if this_peak_group_transition_number >= min_transition:
            peak_group2[rt] = peak_group[rt]
            num_good_transition += 1

    return num_good_transition, peak_group2


def compute_fragment_ranking_expert(this_tg, rt_list, int_list, peak_group, reference_sample_peak_group_info):

    # compute the correlation of the transition intensity of the peak group
    # between this sample and the reference sample

    number_highly_correlated_peak_groups = 0

    peak_group2 = defaultdict(list)

    reference_sample_peak_area = reference_sample_peak_group_info['peak_group_area']

    r_value = {}

    for rt in peak_group:

        this_sample_intensity = compute_transition_intensity(rt, rt_list, int_list)
        x = []
        y = []
        for this_transition in reference_sample_peak_area.keys():
            if this_transition != this_tg:
                x.append(reference_sample_peak_area[this_transition])
                y.append(this_sample_intensity[this_transition])
        r_value[rt] = np.corrcoef(x, y)[0][1]  # NOTE: here is R value, not R **2

    max_r_value = max(r_value.values())

    for rt in peak_group:

        if r_value[rt] >= max_r_value - 0.3 and r_value[rt] >= 0.4:  # have to be >0.4
            peak_group2[rt] = peak_group[rt]
            number_highly_correlated_peak_groups += 1

    return number_highly_correlated_peak_groups, peak_group2


def examine_ms1_expert(this_tg, peak_group):
    # select only the peak groups with MS1
    num_good_ms1 = 0

    peak_group2 = defaultdict(list)

    for rt in peak_group:
        if this_tg in peak_group[rt]:
            peak_group2[rt] = peak_group[rt]
            num_good_ms1 += 1

    return num_good_ms1, peak_group2


def find_best_peak_shape_expert(this_tg, rt_list, int_list, peak_group, reference_sample_peak_group_info):

    peak_group2 = defaultdict(list)

    best_peak_score = -999999.0
    best_peak_rt = 0.0
    reference_sample_peak_width = reference_sample_peak_group_info['peak_group_rt_right'] - \
        reference_sample_peak_group_info['peak_group_rt_left']

    for rt in peak_group:
        peak_intensity = compute_transition_intensity(rt, rt_list, int_list)
        left_intensity, right_intensity = compute_transition_intensity_boundary(
            rt, rt_list, int_list, reference_sample_peak_width)

        peak_shape_score = compute_peak_shape_score(left_intensity, peak_intensity, right_intensity)

        if peak_shape_score > best_peak_score:
            best_peak_score = peak_shape_score
            best_peak_rt = rt

    peak_group2[best_peak_rt] = peak_group[best_peak_rt]

    return peak_group2


def examine_peak_shape_expert(rt, intensity):

    rt_left = rt[0]
    int_left = intensity[0]
    rt_right = rt[-1]
    int_right = intensity[-1]

    rt_apex = -1.0
    int_apex = 0.0
    for rt0, int0 in zip(rt, intensity):
        if int0 > int_apex:
            rt_apex = rt0
            int_apex = int0

    fold_apex_to_left = int_apex / (int_left + 1)
    fold_apex_to_right = int_apex / (int_right + 1)
    fold_left_to_right = (int_left + 1) / (int_right + 1)

    if_good_peak = 0
    if fold_apex_to_left >= 3 and fold_apex_to_right >= 3:  # RULE x
        if_good_peak = 1

    if_interference = 0
    if (fold_left_to_right >= 10 or fold_apex_to_right >= 10) and fold_left_to_right > 10:  # RULE x
        if_interference = 1

    # print "fold_apex_to_left", fold_apex_to_left, "fold right", fold_apex_to_right

    return if_good_peak, if_interference


def get_top_1_transition(rt_list, int_list):
    top_transition = ''
    highest_int = 0
    for this_transition in rt_list.keys():
        sum_intensity = sum(int_list[this_transition])
        if sum_intensity > highest_int:
            highest_int = sum_intensity
            top_transition = this_transition

    return top_transition


def examine_ms2_peak_group_shape_expert(rt_list, int_list):

    # check the top3 transition, if two are good peak, then pass
    top1_transition = get_top_1_transition(rt_list, int_list)

    # print "top 1 transition", top1_transition

    # RULE x. if the peak group is good, depends only on the top 1 transition
    ms2_if_good_peak_group_top1, ms2_if_interference_top1 = examine_peak_shape_expert(
        rt_list[top1_transition], int_list[top1_transition])

    # RULE x. For all the transitions, check if they are interference signal
    ms2_if_good_peak_group = {}
    ms2_if_interference = {}
    for this_transition in rt_list.keys():
        ms2_if_good_peak_group[this_transition], ms2_if_interference[this_transition] = \
            examine_peak_shape_expert(rt_list[this_transition], int_list[this_transition])
        # print ms2_if_interference[this_transition]

    return ms2_if_good_peak_group_top1, ms2_if_interference


def examine_cut_out_peak_group_expert(rt_list, int_list):
    correlation_of_transition_intensity = []
    rt_values = rt_list[rt_list.keys()[0]]
    for i in range(0, len(rt_values)-1, 1):
        rt1_transition_int = []
        rt2_transition_int = []
        for this_transition in rt_list.keys():
            rt1_transition_int.append(int_list[this_transition][i])
            rt2_transition_int.append(int_list[this_transition][i + 1])
        r_value = np.corrcoef(rt1_transition_int, rt2_transition_int)[0][1]
        correlation_of_transition_intensity.append(r_value)
    if_peak_group_cut_out = 0
    if_ms2_peak_group = 0
    # RULE x, very stringent
    if sum(correlation_of_transition_intensity) / len(correlation_of_transition_intensity) > 0.9:
        if_peak_group_cut_out = 1
        if_ms2_peak_group = 1
    return if_peak_group_cut_out, if_ms2_peak_group


def quant_peak_group_expert(display_data):

    # find whether a peak group exists
    # compute peak group area for ms1 and ms2

    ms2_area_background = 99999999.0

    # double check whether peak group exists
    for sample_id in display_data.keys():

        display_data[sample_id]['if_ms1_good_peak'], display_data[sample_id]['if_ms1_interference_peak'] = \
            examine_peak_shape_expert(
                display_data[sample_id]['ms1_rt'], display_data[sample_id]['ms1_int'])

        display_data[sample_id]['ms1_area'] = chrom_peak_boundary.compute_peak_area(display_data[sample_id]['ms1_rt'],
                                                                                    display_data[sample_id][
                                                                                        'ms1_int'],
                                                                                    display_data[sample_id][
                                                                                        'ms1_rt'][0],
                                                                                    display_data[sample_id]['ms1_rt'][-1])

        # ms2 data

        display_data[sample_id]['if_ms2_peak_group'], display_data[sample_id]['if_ms2_interference'] = \
            examine_ms2_peak_group_shape_expert(
                display_data[sample_id]['ms2_rt'], display_data[sample_id]['ms2_int'])

        # RULE x: if ms1 is good peak, sometimes ms2 are not good. In this case,
        # force the ms2 as good peak
        if display_data[sample_id]['if_ms1_good_peak'] == 1:
            display_data[sample_id]['if_ms2_peak_group'] = 1

        # RULE x: if ms2 does not have good peak group, sometimes because the peak
        # is cut out due to chrom limit. rescue them in this case
        if_peak_group_cut_out = 0
        if display_data[sample_id]['if_ms2_peak_group'] == 0:
            if_peak_group_cut_out, display_data[sample_id]['if_ms2_peak_group'] = \
                examine_cut_out_peak_group_expert(
                    display_data[sample_id]['ms2_rt'], display_data[sample_id]['ms2_int'])

        # compute area, remove the interference signal
        for this_transition in display_data[sample_id]['ms2_rt'].keys():
            if display_data[sample_id]['if_ms2_interference'][this_transition] == 0:

                display_data[sample_id]['ms2_area'][this_transition] = \
                    chrom_peak_boundary.compute_peak_area(display_data[sample_id]['ms2_rt'][this_transition],
                                                          display_data[sample_id][
                                                              'ms2_int'][this_transition],
                                                          display_data[sample_id][
                                                              'ms2_rt'][this_transition][0],
                                                          display_data[sample_id]['ms2_rt'][this_transition][-1])
            else:
                # if it is interference trace, assign the area as 0
                display_data[sample_id]['ms2_area'][this_transition] = 0

            if display_data[sample_id]['ms2_area'][this_transition] < ms2_area_background:
                ms2_area_background = display_data[sample_id]['ms2_area'][this_transition]

    return display_data, ms2_area_background


def compute_peak_shape_score(left_intensity, peak_intensity, right_intensity):

    transition_score = {}

    for this_transition in left_intensity.keys():
        score = math.log(
            (peak_intensity[this_transition]+1)/(left_intensity[this_transition]+1), 10)
        score += math.log((peak_intensity[this_transition]+1) /
                          (right_intensity[this_transition]+1), 10)
        transition_score[this_transition] = score

    return transition_score


def reference_sample_peak_number_expert(rt_list, int_list, rt_left, rt_right):
    peaks_rt = peaks.find_peaks_for_a_sample(rt_list, int_list)
    num_transition_with_single_peak = 0
    if_bad_sample = 0
    peaks = []
    for this_transition in peaks_rt:
        for rt in peaks_rt[this_transition]:
            if rt_left < rt < rt_right:
                peaks.append(rt)
    if max(peaks) - min(peaks) > 200:
        if_bad_sample = 1

    return if_bad_sample


def compute_transition_intensity_boundary(rt0, rt_list, int_list, reference_sample_peak_width):
    # based on a rt value, find the best match rt value, then find out the chrom peak intensity
    left_intensity = {}
    right_intensity = {}
    for this_transition in rt_list.keys():

        rt_distance = 9999
        this_intensity = 0.0
        for rt, i in zip(rt_list[this_transition], int_list[this_transition]):
            if abs(rt - rt0 + 0.5 * reference_sample_peak_width) < rt_distance:
                rt_distance = abs(rt - rt0 + 0.5 * reference_sample_peak_width)
                this_intensity = i
        left_intensity[this_transition] = this_intensity

        rt_distance = 9999
        this_intensity = 0.0
        for rt, i in zip(rt_list[this_transition], int_list[this_transition]):
            if abs(rt - rt0 - 0.5 * reference_sample_peak_width) < rt_distance:
                rt_distance = abs(rt - rt0 - 0.5 * reference_sample_peak_width)
                this_intensity = i
        right_intensity[this_transition] = this_intensity

    return left_intensity, right_intensity


def compute_displayed_peak_boundary(left_rt, right_rt, if_peak_group_found, best_sample_id):

    # step 1. get the width of all peaks in all samples.
    # Get the max width so that this is used for the display of all samples
    # 2015.8.3. Here has a bug. In 10_AAAAAAALQAK_2, for guot_L130625_002_SW , there is no good peak group, so the peak boundary is wrong here 140s, while all the others are ~50s.
    # Therefore I need to remove the bad samples for this peak boundary computation, can not simply pick up max value
    # still the following does not work. because in some peptide, maybe most samples have no good peak group.
    # Then I use the range for the best sample.
    # peak_width ={}
    # for sample in left_rt.keys():
    #     if if_peak_group_found[sample] == 1:
    #         peak_width[sample] = right_rt[sample] - left_rt[sample]
    #     else:
    #         print "this sample is not good ", sample
    # max_width = max(peak_width.values())

    peak_width = right_rt[best_sample_id] - left_rt[best_sample_id] + 10

    # print "peak width are ", peak_width

    # step 2. check the peak width of every sample, if the width is lower than the max width,
    # divide the extra width into two parts, and use each to expand on left and right sides
    left_rt2 = {}
    right_rt2 = {}
    for sample in left_rt.keys():
        left_rt2[sample] = 0.5 * (right_rt[sample] + left_rt[sample]) - 0.5 * peak_width
        right_rt2[sample] = 0.5 * (right_rt[sample] + left_rt[sample]) + 0.5 * peak_width

    return left_rt2, right_rt2


# used for obtaining the max intensity for MS1 and MS2 for all transitions for all samples for a tg_id
# the max values are used to normalized the intensities for plotting.
def get_max_intensity_MS1_MS2(in_file, tg_id, sample_names):
    max_intensity_ms1 = 0.0
    max_intensity_ms2 = 0.0

    # get max intensity for MS1 and MS2
    with gzip.open(in_file, 'rb') as IN_FILE:
        reader = csv.DictReader(IN_FILE, delimiter="\t")
        for row in reader:
            for sample_id in sample_names:
                if row['transition_name'] == tg_id:
                    MS1_int_list = map(float, row[sample_id + '_int'].split(","))
                    if (max_intensity_ms1 < max(MS1_int_list)):
                        max_intensity_ms1 = max(MS1_int_list)

                if row['transition_group_id'] == tg_id and row['transition_name'] != tg_id:
                    this_transition_int = map(float, row[sample_id + "_int"].split(","))
                    if (max_intensity_ms2 < max(this_transition_int)):
                        max_intensity_ms2 = max(this_transition_int)

    return max_intensity_ms1, max_intensity_ms2


def multiple_line_text(label):
    label2 = []
    for i in range(len(label)):
        if i % 5 == 0 and i > 0:
            label2.append("\n" + label[i])
        else:
            label2.append(label[i])
    return ''.join(label2)


def write_out_quant_data_expert(quant_data_write_out,
                                display_data,
                                reference_sample_peak_group_info,
                                this_tg, sample_label, peak_boundary_display,
                                ms2_area_background):

    for sample_id in display_data.keys():

        if_reference_sample = 0
        if sample_id == reference_sample_peak_group_info['sample_id']:
            if_reference_sample = 1

        transition_list = [i for i in display_data[sample_id]['ms2_rt'].keys() if i != this_tg]

        transition_list_csv_string = ','.join(transition_list)

        # ms2 area
        transition_ms2_area = []
        for this_transition in transition_list:
            # RULE: if ms2 trace is a interference, do not use for quantification
            if display_data[sample_id]['if_ms2_interference'][this_transition] == 0:
                transition_ms2_area.append(display_data[sample_id]['ms2_area'][this_transition])

        transition_ms2_area_csv_string = ','.join(map(str, transition_ms2_area))
        transition_ms2_area_sum = sum(transition_ms2_area)

        transition_ms2_area_sum2 = transition_ms2_area_sum  # consider whether peak group exists
        if display_data[sample_id]['if_ms2_peak_group'] == 0:
            transition_ms2_area_sum2 = ms2_area_background

        quant_data_write_out.append((this_tg, transition_list_csv_string, sample_id,
                                     sample_label[sample_id], if_reference_sample,
                                     peak_boundary_display['left'][
                                         sample_id], peak_boundary_display['right'][sample_id],
                                     display_data[sample_id][
                                         'if_ms2_peak_group'], transition_ms2_area_csv_string,
                                     transition_ms2_area_sum, transition_ms2_area_sum2, ms2_area_background,
                                     display_data[sample_id]['if_ms1_good_peak'],
                                     display_data[sample_id]['ms1_area']
                                     ))

    return quant_data_write_out
