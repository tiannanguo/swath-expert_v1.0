__author__ = 'Tiannan Guo, ETH Zurich 2015'

import csv
import numpy as np
import operator


def median(lst):
    return np.median(np.array(lst))


def computer_other_sample_area_ratio_median_after_removing_bad_shaped_fragment(i, tg, sample, display_data):
    ratio_list = {}
    i.seek(0)
    r = csv.DictReader(i, delimiter="\t")
    for row in r:
        if row['transition_group_id'] == tg:
            fragment_id = row['fragment_id']
            rt_list = display_data[tg][sample]['ms2']['rt_list'][fragment_id]
            i_list = display_data[tg][sample]['ms2']['i_list'][fragment_id]
            if_good_shape_fragment = check_if_displayed_peak_a_good_one(rt_list, i_list, 1, 1.5)
            if if_good_shape_fragment == 1:
                this_ratio = row[sample + '_ratio_to_ref']
                if this_ratio != 'NA':
                    ratio_list[fragment_id] = float(this_ratio)
    if len(ratio_list.values()) >= 1:
        return median(ratio_list.values())
    else:
        return 'NA'


def compute_peptide_intensity_based_on_median_ratio_of_fragments(quant_file_peptides, quant_file_fragments, sample_id, ref_sample_data, display_data):

    # use information from quant_file_fragments to write out the quant table
    with open(quant_file_peptides, 'wb') as o, open(quant_file_fragments, 'rb') as i:

        w = csv.writer(o, delimiter="\t")

        title = write_title_for_peptide_quant_file(sample_id)

        w.writerow(title)

        for tg in display_data.keys():

            ref_sample_id = ref_sample_data[tg].sample_name

            ref_sample_all_fragments_area = get_ref_sample_all_fragments_area(
                display_data, tg, ref_sample_id)

            data_list = [tg]

            for sample in sample_id:

                # print sample
                if sample != ref_sample_id:

                    if sample == 'nci1':
                        pass

                    this_sample_i_ratio_median = computer_other_sample_area_ratio_median_after_removing_bad_shaped_fragment(
                        i, tg, sample, display_data)

                    this_sample_i = 'NA'

                    if this_sample_i_ratio_median != 'NA':

                        this_sample_i = this_sample_i_ratio_median * ref_sample_all_fragments_area

                    data_list.append(this_sample_i)

                else:
                    data_list.append(ref_sample_all_fragments_area)

            data_list2 = fill_in_background_value(data_list)
            w.writerow(data_list2)
            # w.writerow(data_list)

    return display_data


def get_ref_sample_all_fragments_area(display_data, tg, ref_sample_id):

    i_sum = 0

    for fragment in display_data[tg][ref_sample_id]['ms2']['area'].keys():
        this_i = display_data[tg][ref_sample_id]['ms2']['area'][fragment]
        i_sum += this_i

    return i_sum


def compute_peptide_intensity(display_data, sample_id, ref_sample_data, quant_file_peptides):

    # use information from quant_file_fragments to write out the quant table
    with open(quant_file_peptides, 'wb') as o:

        w = csv.writer(o, delimiter="\t")

        title = write_title_for_peptide_quant_file(sample_id)
        w.writerow(title)

        for tg in display_data.keys():

            ref_sample_id = ref_sample_data[tg].sample_name
            ref_sample_top1_fragment, ref_sample_top1_fragment_i = get_ref_sample_top1_good_shape_fragment_i(
                display_data, tg, ref_sample_id)

            data_list = []
            data_list.append(tg)

            for sample in sample_id:

                # print sample
                if sample != ref_sample_id:

                    # if sample == 'gold36':
                    #     pass

                    other_sample_top1_fragments_i_ratio = compute_other_sample_top1_fragments_i(
                        ref_sample_top1_fragment, display_data, tg, sample, ref_sample_id)
                    other_sample_i = 'NA'

                    if other_sample_top1_fragments_i_ratio != 'NA':

                        other_sample_i = other_sample_top1_fragments_i_ratio * \
                            ref_sample_top1_fragment_i

                    data_list.append(other_sample_i)

                else:
                    data_list.append(ref_sample_top1_fragment_i)

            data_list2 = fill_in_background_value(data_list)
            w.writerow(data_list2)

    return display_data


def fill_in_background_value(data_list):

    i_list = []

    for i in data_list[1:]:
        if i != 'NA':
            i_list.append(float(i))

    if len(i_list) > 1:
        min_i = min(i_list) * 0.1  # arbitrary, can be changed

        data_list2 = []

        data_list2.append(data_list[0])

        for i in data_list[1:]:
            if i == 'NA':
                data_list2.append(min_i)
            else:
                data_list2.append(i)

        return data_list2
    else:
        return data_list


def compute_other_sample_top1_fragments_i(ref_sample_top1_fragment, display_data, tg, sample, ref_sample_id):

    top1_ratio = display_data[tg][sample]['ms2']['ratio_to_ref'][ref_sample_top1_fragment]

    # 2015.12.16. revise. not to check if_good_pg

    # if_good_pg = compute_fragment_correlation(display_data[tg][sample]['ms2']['area'],
    #                                             display_data[tg][ref_sample_id]['ms2']['area'],
    #                                             ref_sample_top1_fragment)
    # if if_good_pg == 1:

    # return top1_ratio
    # else:
    #     return 'NA'

    return top1_ratio


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
        if_good_shape = check_if_displayed_peak_a_good_one(this_rt_list, this_i_list, 1, 1.5)

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
            ref_sample_peptide_i = compute_peak_group_area_for_reference_sample(
                fragment_list, display_data, tg, ref_sample)

            # check each of the fragment
            for fragment in fragment_list:

                data_list, display_data = compute_quant_data_list_for_a_fragment(
                    fragment, tg, sample_id, ref_sample, display_data, ref_sample_peptide_i)

                w.writerow(data_list)

    return display_data


def compute_peak_group_area_for_reference_sample(fragment_list, display_data, tg, ref_sample):

    ref_sample_peptide_i = 0

    for fragment in fragment_list:
        ref_sample_peptide_i += display_data[tg][ref_sample]['ms2']['area'][fragment]

    return ref_sample_peptide_i


def check_if_displayed_peak_a_good_one(rt_list, i_list, if_found_peak, fold_change_threshold):

    # check if the peak is a good one
    if_good = 0

    if len(i_list) > 5:
        point_left_i = i_list[0]
        point_apex_i = max(i_list)  # i_list[len(rt_list) / 2]
        point_right_i = i_list[-1]

        # only check the left and right, do not check the apex because for many weak fragments, the apex is not obvious
        # the left_i and right_i should be both < apex_i, and they are similar
        fold_change_left = point_left_i / (point_apex_i + 0.1)
        fold_change_right = point_right_i / (point_apex_i + 0.1)
        if_good_fold_change_left = check_peak_i_fold_change(fold_change_left, fold_change_threshold)
        if_good_fold_change_right = check_peak_i_fold_change(
            fold_change_right, fold_change_threshold)

        fold_change_left_to_right = (point_left_i + 0.1) / (point_right_i + 0.1)

        # if left intensity and right intensity are comparable, then it is a good
        # shape (note: some weak fragment may be flat)
        if_good = check_left_and_right_only(
            fold_change_left_to_right, if_found_peak, fold_change_threshold)

        # some strong fragment may be cut at one end or both. If the apex is
        # higher than the boundaries, still considered as good shape
        if if_good == 0:
            if_good = check_apex_and_boundary(
                if_good_fold_change_left, if_good_fold_change_right, if_found_peak)

    return if_good


def check_if_displayed_peak_a_good_one_ms1(rt_list, i_list, if_found_peak, fold_change_threshold):

    # ms1 is more strict than ms2
    # check if the peak is a good one
    if_good = 0

    if len(i_list) > 5:
        point_left_i = i_list[0]
        point_apex_i = max(i_list)  # i_list[len(rt_list) / 2]
        point_right_i = i_list[-1]

        # only check the left and right, do not check the apex because for many weak fragments, the apex is not obvious
        # the left_i and right_i should be both < apex_i, and they are similar
        fold_change_left = point_left_i / (point_apex_i + 0.1)
        fold_change_right = point_right_i / (point_apex_i + 0.1)
        if_good_fold_change_left = check_peak_i_fold_change(fold_change_left, fold_change_threshold)
        if_good_fold_change_right = check_peak_i_fold_change(
            fold_change_right, fold_change_threshold)

        fold_change_left_to_right = (point_left_i + 0.1) / (point_right_i + 0.1)

        # if left intensity and right intensity are comparable, then it is a good
        # shape (note: some weak fragment may be flat)
        if_good1 = check_left_and_right_only(
            fold_change_left_to_right, if_found_peak, fold_change_threshold)

        # some strong fragment may be cut at one end or both. If the apex is
        # higher than the boundaries, still considered as good shape
        if_good2 = check_apex_and_boundary(
            if_good_fold_change_left, if_good_fold_change_right, if_found_peak)

        if_good = if_good1 * if_good2

    return if_good


def check_apex_and_boundary(if_good_fold_change_left, if_good_fold_change_right, if_found_peak):

    if if_good_fold_change_left == 1 and if_good_fold_change_right == 1 and if_found_peak == 1:
        return 1
    else:
        return 0


def check_left_and_right_only(fold_change_left_to_right, if_found_peak, fold_change_threshold):

    if 1.0 / float(fold_change_threshold) <= fold_change_left_to_right <= float(fold_change_threshold) and if_found_peak == 1:
        return 1
    else:
        return 0


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

            rt_list = display_data[tg][sample]['ms2']['rt_list'][fragment]
            i_list = display_data[tg][sample]['ms2']['i_list'][fragment]
            if_good_shape_fragment = check_if_displayed_peak_a_good_one(rt_list, i_list, 1, 1.5)

            if if_good_shape_fragment == 1:

                this_area = float(display_data[tg][sample]['ms2']['area'][fragment])
                ref_area = float(display_data[tg][ref_sample]['ms2']['area'][fragment])

                this_ratio = this_area / (ref_area + 0.1)

                display_data[tg][sample]['ms2']['area_refined'][fragment] = this_area
                display_data[tg][sample]['ms2']['ratio_to_ref'][fragment] = this_ratio
                data_list.append(this_area)
                data_list.append(this_ratio)

            else:

                display_data[tg][sample]['ms2']['area_refined'][fragment] = 'NA'
                display_data[tg][sample]['ms2']['ratio_to_ref'][fragment] = "NA"
                data_list.append('NA')
                data_list.append('NA')

        # process reference sample
        else:

            this_area = float(display_data[tg][ref_sample]['ms2']['area'][fragment])
            display_data[tg][ref_sample]['ms2']['area_refined'][fragment] = this_area
            display_data[tg][ref_sample]['ms2']['ratio_to_ref'][fragment] = 1.0

            data_list.append(this_area)
            data_list.append(1.0)

    return data_list, display_data
