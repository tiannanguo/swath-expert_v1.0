__author__ = 'guot'


import numpy as np
from collections import defaultdict
import peaks
import parameters
import data_holder


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
        peak_rt = ref_sample_data[tg].peak_rt
        peak_rt_found = ref_sample_data[tg].peak_rt_found

        # get the peak boundary for the reference sample
        fragments = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments
        i = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_i
        rt_left = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_rt_left
        rt_right = peak_group_candidates[tg][reference_sample][peak_rt_found].matched_fragments_rt_right

        ref_sample_rt_left, ref_sample_rt_right = get_peak_group_boundary(fragments, i, rt_left, rt_right)

        ref_sample_data[tg].read_peak_boundary(ref_sample_rt_left, ref_sample_rt_right)

        # store the chrom to be displayed for the reference sample
        display_data[tg][reference_sample]['rt_left'] = ref_sample_rt_left
        display_data[tg][reference_sample]['rt_right'] = ref_sample_rt_right

        # use the range to narrow MS1 and fragment rt range
        display_data = refine_reference_sample_rt_range(display_data, chrom_data, tg, reference_sample, peak_rt_found)

    return display_data

def refine_reference_sample_rt_range(display_data, chrom_data, tg, reference_sample, peak_rt_found):

    # process MS1 and MS2 together
    for fragment in chrom_data[tg][reference_sample].keys():
        rt_list = chrom_data[tg][reference_sample][fragment].rt_list
        i_list = chrom_data[tg][reference_sample][fragment].i_list
        display_data[tg][reference_sample]['rt_list'][fragment], display_data[tg][reference_sample]['i_list'][fragment] = \
            get_chrom_range(display_data[tg][reference_sample]['rt_left'], display_data[tg][reference_sample]['rt_right'], rt_list, i_list)

    return display_data

def get_peak_group_boundary(fragments, i, rt_left, rt_right):

    # sort by decreasing intensity
    i2 = sorted(i, reverse=1)
    fragments2 = [x for (y, x) in sorted(zip(i, fragments), reverse=1)]
    rt_left2 = [x for (y, x) in sorted(zip(rt_left, fragments), reverse=1)]
    rt_right2 = [x for (y, x) in sorted(zip(rt_right, fragments), reverse=1)]

    rt_left3 = -1
    rt_right3 = -1
    # check the highest fragment first, if a reasonable peak boundary is found, use it. Otherwise, decending the fragment untill find a reasonable boundary
    for fragment, i, rt_left0, rt_right0 in zip(fragments2, i2, rt_left2, rt_right2):
        if rt_right0 - rt_left0 > parameters.MAX_PEAK_WIDTH:
            continue
        else:
            rt_left3 = rt_left0
            rt_right3 = rt_right0
            break
    # if no peak boundary found, use the one from highest peak
    if rt_left3 < 0:
        rt_left3 = rt_left2[0]
        rt_right3 = rt_right2[0]

    return rt_left3, rt_right3

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

def compute_peak_area_for_all(d):
    for tg in d.keys():
        for sample in d[tg]['display_pg'].keys():
            # ms1
            # the displayed rt_list is already cut at the peak boudnary used for display
            area_ms1 = compute_peak_area2(d[tg]['display_pg'][sample]['ms1']['rt_list'], d[tg]['display_pg'][sample]['ms1']['i_list'])
            d[tg]['display_pg'][sample]['area']['ms1'] = area_ms1

            # ms2
            for fragment in d[tg]['display_pg'][sample]['rt_list'].keys():
                area_ms2 = compute_peak_area2(d[tg]['display_pg'][sample]['rt_list'][fragment], d[tg]['display_pg'][sample]['i_list'][fragment])
                d[tg]['display_pg'][sample]['area']['ms2'][fragment] = area_ms2
    return d


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
    rt_list2 = []
    i_list2 = []
    area = 0
    for i in range(1, len(rt_list), 1):
        area += 0.5 * (i_list[i] + i_list[i-1]) * (rt_list[i] - rt_list[i-1])
        rt_list2.append[rt_list[i]]
        i_list2.append[i_list[i]]
    return area



def get_peak_boundary(rt_list, i_list, peak_rt):

    peak_rt_left = rt_list[0]
    peak_rt_right = rt_list[-1]

    max_int = -1.0
    peak_index = -1
    for i in range(len(rt_list)):
        if abs (rt_list[i] - peak_rt ) < 3:
            max_int = i_list[i]
            peak_index = i
            break

    for j in range(peak_index, len(rt_list)):
        #find the boundary, intensity < 10% of max intensity
        if i_list[j] < max_int * 0.1:
            peak_rt_right = rt_list[j]
            break

    for j in range(peak_index, -1, -1):
        #find the boundary, intensity < 10% of max intensity
        if i_list[j] < max_int * 0.1:
            peak_rt_left = rt_list[j]
            break

    #extent 5 seconds
    if peak_rt_left - rt_list[0] > 5:
        peak_rt_left -= 5
    else:
        peak_rt_left = rt_list[0]
    if rt_list[-1] - peak_rt_right > 5:
        peak_rt_right += 5
    else:
        peak_rt_right = rt_list[-1]

    return peak_rt_left, peak_rt_right





