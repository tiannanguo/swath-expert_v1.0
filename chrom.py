__author__ = 'guot'


import numpy as np
from collections import defaultdict
import peaks
import parameters
import data_holder
import peak_groups
import chrom
import swath_quant


# def peak_group_boundary_all_samples (peak_group_info, reference_sample_id):
#
#     reference_sample_peak_width = peak_group_info[reference_sample_id]['peak_group_rt_right'] - \
#         peak_group_info[reference_sample_id]['peak_group_rt_left']
#
#     rt_boundary = defaultdict(dict)
#
#     for sample_id in peak_group_info.keys():
#
#         if 'peak_group_rt_left' in peak_group_info[sample_id].keys():
#
#             rt_boundary['left'][sample_id] = peak_group_info[sample_id]['peak_group_rt_apex'] - \
#                 0.5 * reference_sample_peak_width
#
#             rt_boundary['right'][sample_id] = rt_boundary['left'][sample_id] + \
#                 reference_sample_peak_width
#
#         else:
#             rt_boundary['left'][sample_id] = peak_group_info[reference_sample_id]['peak_group_rt_left']
#             rt_boundary['right'][sample_id] = peak_group_info[reference_sample_id]['peak_group_rt_right']
#
#     return rt_boundary


def compute_reference_sample_peak_boundary(ref_sample_data, chrom_data, peptide_data, peak_group_candidates):

    display_data = data_holder.Nested_dict()

    for tg in chrom_data.keys():

        reference_sample = ref_sample_data[tg].sample_name
        peak_rt_found = ref_sample_data[tg].peak_rt_found

        # get the peak boundary for the reference sample
        fragments = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments
        i = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_i
        rt_left_list = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_peak_rt_left
        rt_right_list = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_peak_rt_right

        ref_sample_rt_left, ref_sample_rt_right = get_peak_group_boundary(fragments, i, rt_left_list, rt_right_list, peak_rt_found)

        # if the ref sample peak is too narrow, < 30 sec, extend 5 sec at both ends
        if ref_sample_rt_right - ref_sample_rt_left < 30:
            ref_sample_rt_left -= 3
            ref_sample_rt_right += 3

        ref_sample_data[tg].read_peak_boundary(ref_sample_rt_left, ref_sample_rt_right)

        # store the chrom to be displayed for the reference sample
        display_data[tg][reference_sample]['rt_left'] = ref_sample_rt_left
        display_data[tg][reference_sample]['rt_right'] = ref_sample_rt_right

        # use the range to narrow MS1 and fragment rt range
        display_data = refine_reference_sample_rt_range(display_data, chrom_data, tg, reference_sample, peak_rt_found)

        display_data, peak_group_candidates, chrom_data = further_refine_ref_sample_fragments(display_data, peak_group_candidates, tg, reference_sample, chrom_data)

    return display_data, peak_group_candidates, chrom_data

def further_refine_ref_sample_fragments(display_data, peak_group_candidates, tg, reference_sample, chrom_data):

    # sometimes, some fragment in the ref sample is too broad
    # it shows a broad peak group in chrom
    # but in the defined boundary, it is not good
    # remove these fragments in all samples
    for fragment in display_data[tg][reference_sample]['ms2']['rt_list'].keys():
        this_rt_list = display_data[tg][reference_sample]['ms2']['rt_list'][fragment]
        this_i_list = display_data[tg][reference_sample]['ms2']['i_list'][fragment]
        if_good_fragment = swath_quant.check_if_displayed_peak_a_good_one(this_rt_list, this_i_list, 1)

        if if_good_fragment == 0:

            # if not good, then remove the fragment
            del display_data[tg][reference_sample]['ms2']['i_list'][fragment]
            del display_data[tg][reference_sample]['ms2']['rt_list'][fragment]
            peak_group_candidates = delete_fragment_from_peak_group_candidate(peak_group_candidates, fragment)
            chrom_data = delete_fragment_from_chrom_data(chrom_data, fragment)
    return display_data, peak_group_candidates, chrom_data

def delete_fragment_from_chrom_data(chrom_data, fragment):

    for tg in chrom_data.keys():
        for sample in chrom_data[tg].keys():
            del chrom_data[tg][sample][fragment]

    return chrom_data

def delete_fragment_from_peak_group_candidate(peak_group_candidates, fragment):

    for tg in peak_group_candidates.keys():
        for sample in peak_group_candidates[tg].keys():
            for rt in peak_group_candidates[tg][sample].keys():
                fragment_index = get_fragment_index(peak_group_candidates[tg][sample][rt].matched_fragments, fragment)
                if fragment_index == -1:
                    pass
                else:
                    del peak_group_candidates[tg][sample][rt].matched_fragments[fragment_index]
                    del peak_group_candidates[tg][sample][rt].matched_fragments_i[fragment_index]
                    del peak_group_candidates[tg][sample][rt].matched_fragments_peak_rt_left[fragment_index]
                    del peak_group_candidates[tg][sample][rt].matched_fragments_peak_rt_right[fragment_index]
                    del peak_group_candidates[tg][sample][rt].matched_fragments_rt[fragment_index]
                    peak_group_candidates[tg][sample][rt].num_matched_fragments -= 1

    return peak_group_candidates

def get_fragment_index(fragments_list, fragment):
    i2 = -1
    for i in range(len(fragments_list)):
        if fragments_list[i] == fragment:
            i2 = i
            break
    return i2

def get_rt_list_in_a_range(rt_list, rt_left, rt_right):

    rt_list2 = []

    for rt in rt_list:
        if float(rt_left) < rt < float(rt_right):
            rt_list2.append(rt)

    return rt_list2

def refine_reference_sample_rt_range(display_data, chrom_data, tg, reference_sample, peak_rt_found):

    # process MS1 and MS2 together
    for fragment in chrom_data[tg][reference_sample].keys():
        rt_list = chrom_data[tg][reference_sample][fragment].rt_list
        i_list = chrom_data[tg][reference_sample][fragment].i_list
        if fragment == tg:
            display_data[tg][reference_sample]['ms1']['rt_list'], display_data[tg][reference_sample]['ms1']['i_list'] = \
            get_chrom_range(display_data[tg][reference_sample]['rt_left'], display_data[tg][reference_sample]['rt_right'], rt_list, i_list)
        else:
            display_data[tg][reference_sample]['ms2']['rt_list'][fragment], display_data[tg][reference_sample]['ms2']['i_list'][fragment] = \
            get_chrom_range(display_data[tg][reference_sample]['rt_left'], display_data[tg][reference_sample]['rt_right'], rt_list, i_list)

    return display_data

def get_peak_group_boundary(fragments, i, rt_left_list, rt_right_list, peak_rt_found):

    pg_rt_left = -1
    pg_rt_right = -1

    # sort fragments by decreasing intensity
    i2 = sorted(i, reverse=1)
    fragments2 = [x for (y, x) in sorted(zip(i, fragments), reverse=1)]
    rt_left_list2 = [y for (y, x) in sorted(zip(rt_left_list, fragments), reverse=1)]
    rt_right_list2 = [y for (y, x) in sorted(zip(rt_right_list, fragments), reverse=1)]

    # check the highest fragment first, if a reasonable peak boundary is found, use it. Otherwise, descending the fragment until find a reasonable boundary
    for fragment, i, rt_left0, rt_right0 in zip(fragments2, i2, rt_left_list2, rt_right_list2):

        distance_left = peak_rt_found - rt_left0
        distance_right = rt_right0 - peak_rt_found

        if_good_peak_boundary = 0
        if distance_left > 0 and distance_right > 0:
            ratio = float(distance_left) / float(distance_right)
            if parameters.PEAK_BOUNDARY_RT_LEFT_RIGHT_RATIO_TOLERANCE < ratio < 1 / parameters.PEAK_BOUNDARY_RT_LEFT_RIGHT_RATIO_TOLERANCE:
                if_good_peak_boundary = 1

        if rt_right0 - rt_left0 > parameters.MAX_PEAK_WIDTH or if_good_peak_boundary == 0:

            continue

        else:

            peak_half_width = min(distance_left, distance_right)
            pg_rt_left = peak_rt_found - peak_half_width
            pg_rt_right = peak_rt_found + peak_half_width
            break

    # if no peak boundary found, use the one from highest peak
    if pg_rt_left < 0:
        pg_rt_left = rt_left_list2[0]
        pg_rt_right = rt_right_list2[0]
    else:
        pg_rt_left += 3
        pg_rt_right += 3

    return pg_rt_left, pg_rt_right

def get_chrom_range(rt_left, rt_right, rt_list, i_list):
    rt_list2 = []
    i_list2 = []
    for rt, i in zip(rt_list, i_list):
        if rt_left <= rt <= rt_right:
            rt_list2.append(rt)
            i_list2.append(i)
    return rt_list2, i_list2

def find_closest_rt(rt, rt_list):
    rt2 = rt_list[0]
    dif = abs(rt2 - rt)
    for rt0 in rt_list:
        if abs(rt0 - rt) < dif:
            dif = abs(rt0 - rt)
            rt2 = rt0
    return rt2

def compute_peak_area_for_all(display_data):

    for tg in display_data.keys():

        for sample in display_data[tg].keys():
            # ms1
            # the displayed rt_list is already cut at the peak boundary used for display
            area_ms1 = compute_peak_area2(display_data[tg][sample]['ms1']['rt_list'], display_data[tg][sample]['ms1']['i_list'])
            display_data[tg][sample]['ms1']['area'] = area_ms1

            # ms2
            for fragment in display_data[tg][sample]['ms2']['rt_list'].keys():

                area_ms2 = compute_peak_area2(display_data[tg][sample]['ms2']['rt_list'][fragment],
                                              display_data[tg][sample]['ms2']['i_list'][fragment])

                display_data[tg][sample]['ms2']['area'][fragment] = area_ms2

    return display_data


def compute_peak_area(rt_list, i_list, rt_left, rt_right):
    rt_list2 = []
    i_list2 = []
    area = 0
    for i in range(1, len(rt_list), 1):
        if rt_left < rt_list[i] < rt_right :
            area += 0.5 * (i_list[i] + i_list[i-1]) * (rt_list[i] - rt_list[i-1])
            rt_list2.append[rt_list[i]]
            i_list2.append[i_list[i]]
    return rt_list2, i_list2, area

def compute_peak_area2(rt_list, i_list):

    # the rt_list and i_list is already refined. directly compute the area

    area = 0

    for i in range(1, len(rt_list), 1):
        area += 0.5 * (i_list[i] + i_list[i-1]) * (rt_list[i] - rt_list[i-1])

    return area



def get_peak_boundary(rt_list, i_list, peak_rt):

    # find the apex
    peak_index = -1
    for i in range(len(rt_list)):
        if abs(rt_list[i] - peak_rt) < 3:
            peak_index = i
            break

    # find right boundary
    start_index_right = peak_index
    end_index_right = len(rt_list) - 1
    peak_rt_right = get_peak_boundary_worker_right(start_index_right, end_index_right, rt_list, i_list)

    # find left boundary
    start_index_left = peak_index
    end_index_left = 0
    peak_rt_left = get_peak_boundary_worker_left(start_index_left, end_index_left, rt_list, i_list)

    # extent 3 seconds
    if peak_rt_left - rt_list[0] > 3:
        peak_rt_left -= 3
    else:
        peak_rt_left = rt_list[0]

    if rt_list[-1] - peak_rt_right > 3:
        peak_rt_right += 3
    else:
        peak_rt_right = rt_list[-1]

    return peak_rt_left, peak_rt_right



def get_peak_boundary_worker_right(start_index, end_index, rt_list, i_list):

    # find the turning point by checking slope value
    boundary_index = -1

    peak_i = i_list[start_index]

    for j in range(start_index, end_index):
        point_1_i = i_list[j]
        point_2_i = i_list[j + 1]

        i_ratio = point_1_i / (peak_i + 0.001)

        if i_ratio < 0.3: # empirical value
            if point_1_i <= point_2_i or i_ratio < 0.1:
                # this is a turning point
                boundary_index = j
                break

    if boundary_index == -1:  # if no double-confirmed turning point is found, then use the backup
        # if no turning point found, then use the end point
        boundary_index = end_index

    return rt_list[boundary_index]

def get_peak_boundary_worker_left(start_index, end_index, rt_list, i_list):

    # find the turning point by checking slope value
    boundary_index = -1

    peak_i = i_list[start_index]

    for j in range(start_index, end_index, -1):  #end point + 1 because we consider two points here
        point_1_i = i_list[j]
        point_2_i = i_list[j - 1]

        i_ratio = point_1_i / (peak_i + 0.001)

        if i_ratio < 0.3: # empirical value
            if point_1_i <= point_2_i or i_ratio < 0.1:
                # this is a turning point
                boundary_index = j
                break

    if boundary_index == -1:  # if no double-confirmed turning point is found, then use the backup
        # if no turning point found, then use the end point
        boundary_index = end_index

    return rt_list[boundary_index]
