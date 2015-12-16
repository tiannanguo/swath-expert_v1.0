__author__ = 'Tiannan Guo, ETH Zurich 2015'

import csv
import numpy as np
from scipy import stats
import operator


def compute_peptide_intensity(display_data, sample_id, ref_sample_data, quant_file_peptides):

    # use information from quant_file_fragments to write out the quant table
    with open(quant_file_peptides, 'wb') as o:

        w = csv.writer(o, delimiter="\t")

        title = write_title_for_peptide_quant_file(sample_id)
        w.writerow(title)

        for tg in display_data.keys():

            ref_sample_id = ref_sample_data[tg].sample_name
            ref_sample_top1_fragment, ref_sample_top1_fragment_i = get_ref_sample_top1_good_shape_fragment_i(display_data, tg, ref_sample_id)

            data_list = []
            data_list.append(tg)

            for sample in sample_id:

                # print sample
                if sample != ref_sample_id:

                    # if sample == 'gold36':
                    #     pass

                    other_sample_top1_fragments_i_ratio = compute_other_sample_top1_fragments_i(ref_sample_top1_fragment, display_data, tg, sample, ref_sample_id)
                    other_sample_i = 'NA'

                    if other_sample_top1_fragments_i_ratio != 'NA':

                        other_sample_i = other_sample_top1_fragments_i_ratio * ref_sample_top1_fragment_i

                    data_list.append(other_sample_i)

                else:
                    data_list.append(ref_sample_top1_fragment_i)

            data_list2 = fill_in_background_value(data_list)
            w.writerow(data_list2)

    return display_data

def fill_in_background_value(data_list):

    i_list = []

    for i in data_list[1:]:
        if type(i) == float:
            i_list.append(i)

    if len(i_list) > 1:
        min_i = min(i_list) * 0.1 # arbitrary, can be changed

        data_list2 = []

        data_list2.append(data_list[0])

        for i in data_list[1:]:
            if type(i) == float:
                data_list2.append(i)
            else:
                data_list2.append(min_i)

        return data_list2
    else:
        return data_list

def compute_other_sample_top1_fragments_i(ref_sample_top1_fragment, display_data, tg, sample, ref_sample_id):

    top1_ratio = display_data[tg][sample]['ms2']['ratio_to_ref'][ref_sample_top1_fragment]

    if_good_pg = compute_fragment_correlation(display_data[tg][sample]['ms2']['area'],
                                                display_data[tg][ref_sample_id]['ms2']['area'],
                                                ref_sample_top1_fragment)
    if if_good_pg == 1:

        return top1_ratio
    else:
        return 'NA'

def find_top_n_fragment_based_on_area(option, area):

    # sort fragment based on i, and then select the top n fragment

    area2 = sorted(area.items(), key=operator.itemgetter(1), reverse=True)

    if len(area2) >= int(option):
        return area2[int(option) - 1][0]
    else:
        return 'NA'

def find_top_3_fragments_based_on_area(area):

    top_fragments = {}
    for i in range(1, 4):
        top_i = find_top_n_fragment_based_on_area(i, area)
        if top_i != 'NA':
            top_fragments[i] = top_i

    return top_fragments

def compute_fragment_correlation(area1, area2, top1_fragment):

    # area1 is the sample, area2 is reference,
    # it's a good pg if the following two scenarios happen:
    # 1, the top1 fragment in ref is the top1 in the sample
    # 2, else, the top2 and top3 in ref sample is top 3 in the sample

    top_fragments = find_top_3_fragments_based_on_area(area1)

    ref_top_fragments = find_top_3_fragments_based_on_area(area2)

    if_good_pg = 0

    if len(top_fragments) < 3 or len(ref_top_fragments) < 3:
        return if_good_pg
    else:
        if top_fragments[1] == ref_top_fragments[1]:
            if_good_pg = 1
        elif ref_top_fragments[2] in top_fragments.values() or ref_top_fragments[3] in top_fragments.values():
            if_good_pg = 1

    return if_good_pg

def get_ref_sample_top1_good_shape_fragment_i(display_data, tg, ref_sample_id):

    top1_fragment = ''
    top1_i = -1

    for fragment in display_data[tg][ref_sample_id]['ms2']['area'].keys():
        this_i = display_data[tg][ref_sample_id]['ms2']['area'][fragment]
        this_rt_list = display_data[tg][ref_sample_id]['ms2']['rt_list'][fragment]
        this_i_list = display_data[tg][ref_sample_id]['ms2']['i_list'][fragment]
        if_good_shape = check_if_displayed_peak_a_good_one(this_rt_list, this_i_list, 1, 2)

        if this_i > top1_i and if_good_shape == 1:
            top1_fragment = fragment
            top1_i = this_i

    # maybe no peak shape in the reference sample is good, then pick the top fragment
    if top1_i == -1:
        for fragment in display_data[tg][ref_sample_id]['ms2']['area'].keys():
            this_i = display_data[tg][ref_sample_id]['ms2']['area'][fragment]

            if this_i > top1_i:
                top1_fragment = fragment
                top1_i = this_i

    return top1_fragment, top1_i

def write_title_for_fragment_quant_file(sample_id):

    title = ['transition_group_id', 'fragment_id']

    for sample in sample_id:
        title.append(sample + '_area')
        title.append(sample + '_ratio_to_ref')

    return title

def write_title_for_peptide_quant_file(sample_id):

    title = ['transition_group_id']

    for sample in sample_id:
        title.append(sample + '_area')

    return title


def compute_peak_area_for_refined_fragment(display_data, sample_id, ref_sample_data, quant_file_fragments):

    with open(quant_file_fragments, 'wb') as o:

        w = csv.writer(o, delimiter="\t")

        title = write_title_for_fragment_quant_file(sample_id)
        w.writerow(title)

        for tg in display_data.keys():

            ref_sample = ref_sample_data[tg].sample_name

            fragment_list = display_data[tg][ref_sample]['ms2']['area'].keys()

            # compute peak group area of reference sample
            ref_sample_peptide_i = 0
            for fragment in fragment_list:
                ref_sample_peptide_i += display_data[tg][ref_sample]['ms2']['area'][fragment]

            # check each of the fragment
            for fragment in fragment_list:

                data_list, display_data = compute_quant_data_list_for_a_fragment(fragment, tg, sample_id, ref_sample, display_data, ref_sample_peptide_i)

                w.writerow(data_list)

    return display_data


def check_if_displayed_peak_a_good_one(rt_list, i_list, if_found_peak, fold_change_threshold):

    # check if the peak is a good one
    if_good = 0

    if len(i_list) > 5:
        point_left_i = i_list[0]
        point_apex_i = max(i_list) #i_list[len(rt_list) / 2]
        point_right_i = i_list[-1]

        # the left_i and right_i should be both < apex_i, and they are similar
        fold_change_left = point_left_i / (point_apex_i + 0.1)
        fold_change_right = point_right_i / (point_apex_i + 0.1)
        if_good_fold_change_left = check_peak_i_fold_change(fold_change_left, fold_change_threshold)
        if_good_fold_change_right = check_peak_i_fold_change(fold_change_right, fold_change_threshold)

        if if_good_fold_change_left == 1 and if_good_fold_change_right == 1 and if_found_peak == 1:
            if_good = 1

    return if_good

def check_peak_i_fold_change(fold_change, threshold):

    if_good = 0

    if fold_change <= 1.0 / threshold:
        if_good = 1

    return if_good

def compute_quant_data_list_for_a_fragment(fragment, tg, sample_id, ref_sample, display_data, ref_sample_peptide_i):

    data_list = []

    data_list.append(tg)
    data_list.append(fragment)

    for sample in sample_id:

        # other samples
        if sample != ref_sample:

            # if sample == 'gold4':
            #     pass

            this_area = float(display_data[tg][sample]['ms2']['area'][fragment])
            ref_area = float(display_data[tg][ref_sample]['ms2']['area'][fragment])

            this_ratio = this_area / (ref_area + 0.1)

            display_data[tg][sample]['ms2']['area_refined'][fragment] = this_area
            display_data[tg][sample]['ms2']['ratio_to_ref'][fragment] = this_ratio
            data_list.append(this_area)
            data_list.append(this_ratio)

        # process reference sample
        else:

            this_area = float(display_data[tg][ref_sample]['ms2']['area'][fragment])
            display_data[tg][ref_sample]['ms2']['area_refined'][fragment] = this_area
            display_data[tg][ref_sample]['ms2']['ratio_to_ref'][fragment] = 1.0

            data_list.append(this_area)
            data_list.append(1.0)

    return data_list, display_data
