__author__ = 'Tiannan Guo, ETH Zurich 2015'

import csv
import numpy as np


def compute_peptide_intensity(display_data, sample_id, ref_sample_data, quant_file_peptides):

    # use information from quant_file_fragments to write out the quant table
    with open(quant_file_peptides, 'wb') as o:

        w = csv.writer(o, delimiter="\t")

        title = write_title_for_peptide_quant_file(sample_id)
        w.writerow(title)

        for tg in display_data.keys():

            ref_sample_id = ref_sample_data[tg].sample_name
            ref_sample_top1_fragment, ref_sample_top1_fragment_i = get_ref_sample_top1_fragment_i(display_data, tg, ref_sample_id)

            data_list = []
            data_list.append(tg)

            for sample in sample_id:

                # print sample
                if sample != ref_sample_id:

                    if sample == 'gold1':
                        pass

                    other_sample_top1_fragments_i_ratio = compute_other_sample_top1_fragments_i(display_data, tg, sample, ref_sample_id)
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

    min_i = min(i_list) * 0.1 # arbitrary, can be changed

    data_list2 = []

    data_list2.append(data_list[0])
    for i in data_list[1:]:
        if type(i) == float:
            data_list2.append(i)
        else:
            data_list2.append(min_i)

    return data_list2

def compute_other_sample_top1_fragments_i(display_data, tg, sample, ref_sample_id):

    top1_fragment = ''
    top1_i = -1
    top1_ratio = -1

    for fragment in display_data[tg][sample]['ms2']['area'].keys():
        if display_data[tg][sample]['ms2']['ratio_to_ref'][fragment] != 'NA':
            this_i = display_data[tg][sample]['ms2']['area_refined'][fragment]
            if this_i > top1_i:
                top1_fragment = fragment
                top1_i = this_i
                top1_ratio = display_data[tg][sample]['ms2']['ratio_to_ref'][fragment]

    if top1_i < 0:
        top1_ratio = 'NA'

    fragment_cor = compute_fragment_correlation(display_data[tg][sample]['ms2']['area'], display_data[tg][ref_sample_id]['ms2']['area'])

    if fragment_cor < 0.5: # empirical value
        top1_ratio = 'NA'

    return top1_ratio

def compute_fragment_correlation(area1, area2):

    x = []
    y = []

    for fragment in area1.keys():
        x.append(area1[fragment])
        y.append(area2[fragment])

    corr = np.corrcoef(x, y)[0][1]

    return corr



def get_ref_sample_top1_fragment_i(display_data, tg, ref_sample_id):

    top1_fragment = ''
    top1_i = -1

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


def check_if_displayed_peak_a_good_one(rt_list, i_list, if_found_peak):

    # check if the peak is a good one
    if_good = 0

    point_left_rt = rt_list[0]
    point_apex_rt = rt_list[len(rt_list) / 2]
    point_right_rt = rt_list[-1]

    point_left_i = i_list[0]
    point_apex_i = i_list[len(rt_list) / 2]
    point_right_i = i_list[-1]

    max_i = max(i_list)

    # the apex_i should be close to the max intensity
    if_good_apex_i = 0
    if point_apex_i > 0.9 * max_i:
        if_good_apex_i = 1

    # the left_i and right_i should be both < apex_i, and they are similar
    fold_change_left = point_left_i / (point_apex_i + 0.1)
    fold_change_right = point_right_i / (point_apex_i + 0.1)
    if_good_fold_change_left = check_peak_i_fold_change(fold_change_left, 1.2)
    if_good_fold_change_right = check_peak_i_fold_change(fold_change_right, 1.2)

    if if_good_apex_i == 1 and if_good_fold_change_left == 1 and if_good_fold_change_right == 1 and if_found_peak == 1:
        if_good = 1

    return if_good

def check_peak_i_fold_change(fold_change, threshold):

    if_good = 0

    if fold_change <= threshold:
        if_good = 1

    return if_good

def compute_quant_data_list_for_a_fragment(fragment, tg, sample_id, ref_sample, display_data, ref_sample_peptide_i):

    data_list = []

    data_list.append(tg)
    data_list.append(fragment)

    for sample in sample_id:

        # other samples
        if sample != ref_sample:

            if sample == 'gold4':
                pass

            num_good_fragments = 0.0

            # if found a good peak
            if_found_peak = display_data[tg][sample]['ms2']['if_found_peak'][fragment]
            rt_list = display_data[tg][sample]['ms2']['rt_list'][fragment]
            i_list = display_data[tg][sample]['ms2']['i_list'][fragment]
            if_it_is_a_good_peak_shape = check_if_displayed_peak_a_good_one(rt_list, i_list, if_found_peak)

            # if found a good peak then use it for quantification
            if if_it_is_a_good_peak_shape == 1:

                num_good_fragments += 1.0

                this_area = float(display_data[tg][sample]['ms2']['area'][fragment])
                ref_area = float(display_data[tg][ref_sample]['ms2']['area'][fragment])

                this_ratio = this_area / (ref_area + 0.1)

            else:

            # otherwise, write NA value
                this_area = 'NA'
                this_ratio = 'NA'

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
