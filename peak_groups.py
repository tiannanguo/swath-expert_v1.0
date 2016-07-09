__author__ = 'guot'

import chrom
import numpy as np
import parameters
import data_holder
import operator
import swath_quant
import math


def check_if_ms1_peak(chrom_data, tg, sample, rt):

    if_ms1_peak = 0

    if hasattr(chrom_data[tg][sample][tg], 'peak_apex_rt_list'):
        peak_apex_rt_list = chrom_data[tg][sample][tg].peak_apex_rt_list
        peak_apex_i_list = chrom_data[tg][sample][tg].peak_apex_i_list
        if peak_apex_rt_list != 'NA':
            for rt0, i in zip(peak_apex_rt_list, peak_apex_i_list):
                if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:  # 10 s for 2hr gradient
                    if_ms1_peak = 1
                    break

    return if_ms1_peak


def find_matched_fragments(chrom_data, tg, sample, rt):

    # fold_change = 0.1 by default

    matched_fragments = []
    matched_fragments_rt_list = []
    matched_fragments_i_list = []
    matched_fragments_peak_rt_left = []
    matched_fragments_peak_rt_right = []

    for fragment in chrom_data[tg][sample].keys():

        if fragment != tg:

            if hasattr(chrom_data[tg][sample][fragment], 'peak_apex_rt_list'):
                peak_apex_rt_list = chrom_data[tg][sample][fragment].peak_apex_rt_list
                peak_apex_i_list = chrom_data[tg][sample][fragment].peak_apex_i_list

                if peak_apex_rt_list != 'NA':

                    for rt0, i in zip(peak_apex_rt_list, peak_apex_i_list):

                        if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:  # 10 s for 2hr gradient

                            matched_fragments.append(fragment)
                            matched_fragments_rt_list.append(rt0)
                            matched_fragments_i_list.append(i)
                            rt_list = chrom_data[tg][sample][fragment].rt_list
                            i_list = chrom_data[tg][sample][fragment].i_list
                            rt_left, rt_right = chrom.get_peak_boundary(rt_list, i_list, rt0)
                            matched_fragments_peak_rt_left.append(rt_left)
                            matched_fragments_peak_rt_right.append(rt_right)

                            break

    return matched_fragments, matched_fragments_rt_list, matched_fragments_i_list, matched_fragments_peak_rt_left, matched_fragments_peak_rt_right


def find_best_peak_group_based_on_reference_sample(display_data, ref_sample_data, chrom_data, peptide_data, peak_group_candidates, sample_id):

    for tg in chrom_data.keys():

        if tg in ref_sample_data.keys():

            # build a NestedDict object for the peak group from reference sample
            display_data, ref_pg = build_reference_peak_group(
                display_data, ref_sample_data, chrom_data, tg)

            if ref_pg == "NA":
                continue
                print "remvoe tg ", tg, "in find_best_peak_group_based_on_reference_sample"

            for sample in sample_id:

                # print 'sample is ', sample

                if sample == 'BR1':
                    pass

                if sample != ref_sample_data[tg].sample_name:
                    # for each peak group, create a data structure containing all the information above

                    pg = build_other_sample_peak_group(
                        chrom_data, tg, ref_pg, peak_group_candidates, sample)

                    if len(pg) == 0:
                        print 'tg is ', tg, ', sample is ', sample, 'no pg found'

                    pg_best = find_best_match_pg(pg, ref_pg, sample)

                    # BUG:when there is no pg, fill with the best "peak group"!!!!! -? 151202

                    if pg_best == 0:
                        # TODO no peak found. cannot be!!!
                        # use looser criteria to find peak groups and then select the best one
                        # debug debug debug......
                        print 'WARNING:no best peak group!!'

                    else:

                        # write into display_pg
                        display_data[tg][sample]['rt_left'] = pg_best['rt_left']
                        display_data[tg][sample]['rt_right'] = pg_best['rt_right']

                        for fragment in pg_best['ms2']['rt_list'].keys():
                            display_data[tg][sample]['ms2']['rt_list'][
                                fragment] = pg_best['ms2']['rt_list'][fragment]
                            display_data[tg][sample]['ms2']['i_list'][
                                fragment] = pg_best['ms2']['i_list'][fragment]
                            display_data[tg][sample]['ms2']['peak_apex_i'][
                                fragment] = pg_best['ms2']['peak_apex_i'][fragment]
                            display_data[tg][sample]['ms2']['if_found_peak'][
                                fragment] = pg_best['ms2']['if_found_peak'][fragment]

                        display_data[tg][sample]['ms1']['rt_list'] = pg_best['ms1']['rt_list']
                        display_data[tg][sample]['ms1']['i_list'] = pg_best['ms1']['i_list']
                        display_data[tg][sample]['ms1']['peak_apex_i'] = pg_best['ms1']['peak_apex_i']
                        display_data[tg][sample]['ms1']['if_found_peak'] = pg_best['ms1']['if_found_peak']

    return display_data


def get_fragment_intensity_for_peak_group(peaks_rt, peaks_i, rt):

    rt_dif = 100
    i = 0
    if_found_peak = 0

    for rt0, i0, in zip(peaks_rt, peaks_i):
        rt_dif0 = abs(rt0 - rt)
        if rt_dif0 < rt_dif:
            i = i0
            rt_dif = rt_dif0

    if rt_dif < parameters.MAX_RT_TOLERANCE:
        if_found_peak = 1

    return i, if_found_peak


def get_peak_width_left_and_right(pg, rt, ref_pg_width):

    width_left = {}
    width_right = {}
    ratio_right_to_left = {}

    for fragment in pg[rt]['ms2']['rt_left'].keys():
        if pg[rt]['ms2']['rt_left'][fragment] != -1 and pg[rt]['ms2']['rt_right'][fragment] != -1:
            width_left[fragment] = rt - pg[rt]['ms2']['rt_left'][fragment]
            width_right[fragment] = pg[rt]['ms2']['rt_right'][fragment] - rt
            ratio_right_to_left[fragment] = (
                width_right[fragment] + 0.01) / (width_left[fragment] + 0.01)
        else:
            width_left[fragment] = rt - 0.5 * ref_pg_width
            width_right[fragment] = rt + 0.5 * ref_pg_width
            ratio_right_to_left[fragment] = 1.0

    return width_left, width_right, ratio_right_to_left


def reject_outliers(data, m=2):
    data_mean = np.mean(data)
    data_sd = np.std(data)
    data2 = [x for x in data if abs(x - data_mean) < m * data_sd]
    return data2


def get_boundary_for_pg_best(pg, rt, ref_pg_width):

    # in most cases, the best pg is symmetric, then it's simple, use the center +/- half peak width
    # sometimes, the best pg is not symmetric, then proportionally distribute the peak width

    # check the distribution of peak width
    peak_width_on_the_left, peak_width_on_the_right, peak_width_ratio_right_to_left = \
        get_peak_width_left_and_right(pg, rt, ref_pg_width)

    peak_width_ratio_right_to_left2 = reject_outliers(peak_width_ratio_right_to_left.values())

    peak_width_ratio_right_to_left_mean = 1.0

    if len(peak_width_ratio_right_to_left2) > 1:
        peak_width_ratio_right_to_left_mean = np.mean(peak_width_ratio_right_to_left2)

    pg_left = rt - ref_pg_width / (1.0 + peak_width_ratio_right_to_left_mean)
    pg_right = pg_left + ref_pg_width

    return pg_left, pg_right


def get_peak_group_values(pg, rt, ref_pg):

    pg_best = data_holder.NestedDict()

    pg_best['peak_rt'] = rt

    # use peak width from ref_pg to compute the boundary of this peak group
    ref_pg_width = ref_pg['rt_right'] - ref_pg['rt_left']
    pg_best['rt_left'], pg_best['rt_right'] = get_boundary_for_pg_best(pg, rt, ref_pg_width)

    for fragment in pg[rt]['ms2']['rt_list'].keys():

        # get the rt and i list only for the peak group
        pg_best['ms2']['rt_list'][fragment], pg_best['ms2']['i_list'][fragment] = \
            chrom.get_chrom_range(pg_best['rt_left'], pg_best['rt_right'],
                                  pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])

        pg_best['ms2']['peak_apex_i'][fragment], pg_best['ms2']['if_found_peak'][fragment] = get_fragment_intensity_for_peak_group(
            pg[rt]['ms2']['peak_apex_rt_list'][fragment], pg[rt]['ms2']['peak_apex_i_list'][fragment], rt)

    # get the rt and i list only for the peak group
    #sometimes for MS1, the peak_apex_i_list=NA, peak_apex_rt_list=NA, i_list=NA, rt_list=NA
    if pg[rt]['ms1']['rt_list'] == "NA" or pg[rt]['ms1']['i_list'] == "NA":

        pg_best['ms1']['rt_list'] = "NA"
        pg_best['ms1']['i_list'] = "NA"
        pg_best['ms1']['peak_apex_i'] = "NA"
        pg_best['ms1']['if_found_peak'] = "NA"

    else:
        pg_best['ms1']['rt_list'], pg_best['ms1']['i_list'] = \
            chrom.get_chrom_range(
                pg_best['rt_left'], pg_best['rt_right'], pg[rt]['ms1']['rt_list'], pg[rt]['ms1']['i_list'])

        pg_best['ms1']['peak_apex_i'], pg_best['ms1']['if_found_peak'] = get_fragment_intensity_for_peak_group(
            pg[rt]['ms1']['peak_apex_rt_list'], pg[rt]['ms1']['peak_apex_i_list'], rt)

    return pg_best


def find_top_n_fragment(option, ref_pg):

    # sort fragment based on i, and then select the top n fragment

    peak_apex_i = sorted(
        ref_pg['ms2']['peak_apex_i'].items(), key=operator.itemgetter(1), reverse=True)

    if len(peak_apex_i) >= int(option):
        return peak_apex_i[int(option) - 1][0]
    else:
        return 'NA'


def filter_peak_group_top_fragment(n, pg, ref_pg, pg_filtered_rt):

    pg_filtered_rt2 = []

    top_fragment = find_top_n_fragment(n, ref_pg)

    if top_fragment == 'NA':
        return pg_filtered_rt

    for rt in pg_filtered_rt:

        if_peak_found = 0

        if_peak_is_a_top_peak = check_if_a_peak_is_a_top_peak(top_fragment, rt, n, pg)

        if if_peak_is_a_top_peak == 1:
            if_peak_found = 1

        else:
            # sometimes the top1 peak in reference sample is not the best signal (still openswath issue)
            # then check the top2 peak, if top2 peak is ok, still keep it.

            next_fragment = find_top_n_fragment(n + 1, ref_pg)

            if next_fragment != 'NA':

                if_next_peak_is_a_top_peak = check_if_a_peak_is_a_top_peak(next_fragment, rt, n, pg)

                if if_next_peak_is_a_top_peak == 1:
                    if_peak_found = 1

        if if_peak_found == 1:
            pg_filtered_rt2.append(rt)

    return pg_filtered_rt2


def check_if_a_peak_is_a_top_peak(fragment, rt, n, pg):

    if_good = 0

    for rt0 in pg[rt]['ms2']['peak_apex_rt_list'][fragment]:

        if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:

            rank_num = find_fragment_rank_in_a_pg(fragment, rt0, pg)

            # if the fragment rank is high enough, otherwise it's a wrong peak group
            if rank_num - n <= 3:
                if_good = 1

    return if_good


def find_fragment_rank_in_a_pg(fragment, rt, pg):

    rank_num = 1
    fragment_i = pg[rt]['ms2']['peak_apex_i'][fragment]
    # sometimes fragment_i is empty, assign a low value
    if isinstance(fragment_i, float) == 0:
        fragment_i = 0.1
    for this_fragment in pg[rt]['ms2']['peak_apex_i'].keys():
        if this_fragment != fragment:
            # times 1.2 because sometimes the difference is not big. in that case,
            # rank as high as possible
            if pg[rt]['ms2']['peak_apex_i'][this_fragment] > fragment_i * 1.2:
                rank_num += 1

    return rank_num


def filter_peak_group_ms1(pg, pg_filtered_rt):

    pg_filtered_rt2 = []

    for rt in pg_filtered_rt:
        if_peak_found = 0
        for rt0 in pg[rt]['ms1']['peak_apex_rt_list']:
            if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:
                if_peak_found = 1
                break

        if if_peak_found == 1:
            pg_filtered_rt2.append(rt)

    return pg_filtered_rt2


def filter_peak_group_peak_shape(n, pg, ref_pg, pg_filtered_rt):

    pg_filtered_rt2 = []
    top_fragment = find_top_n_fragment(n, ref_pg)

    if top_fragment == 'NA':
        return pg_filtered_rt

    for rt in pg_filtered_rt:

        if_top_fragment_good = check_peak_group_peak_shape(top_fragment, rt, pg)
        if if_top_fragment_good == 1:
            pg_filtered_rt2.append(rt)  # good peak boundary
        else:
            next_fragment = find_top_n_fragment(n + 1, ref_pg)
            if next_fragment != 'NA':
                if_next_fragment_good = check_peak_group_peak_shape(next_fragment, rt, pg)
                if if_next_fragment_good == 1:
                    pg_filtered_rt2.append(rt)  # good peak boundary

    return pg_filtered_rt2


def check_peak_group_peak_shape(fragment, rt, pg):

    if_good = 0

    if "peak_apex_i" in pg[rt]['ms2'].keys():
        if fragment in pg[rt]['ms2']['peak_apex_i'].keys():
            if pg[rt]['ms2']['peak_apex_i'][fragment] > 0:

                peak_intensity_apex = pg[rt]['ms2']['peak_apex_i'][fragment]

                peak_intensity_left = get_intensity_for_closest_rt(
                    pg[rt]['ms2']['rt_left'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])

                peak_intensity_right = get_intensity_for_closest_rt(
                    pg[rt]['ms2']['rt_right'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])

                fold_change_left = float(peak_intensity_apex) / (float(peak_intensity_left) + 1)
                fold_change_right = float(peak_intensity_apex) / (float(peak_intensity_right) + 1)

                if fold_change_left >= parameters.PEAK_SHAPE_FOLD_VARIATION and fold_change_right >= parameters.PEAK_SHAPE_FOLD_VARIATION:
                    if_good = 1

    return if_good


def filter_pg_based_on_peak_boundary(pg_filtered_rt, pg):

    pg_filtered_rt2 = {}

    for rt in pg_filtered_rt:
        num_bad_boundary = 0
        for fragment in pg[rt]['ms2']['rt_left'].keys():
            if pg[rt]['ms2']['rt_left'][fragment] == -1:
                num_bad_boundary += 1
        pg_filtered_rt2[rt] = num_bad_boundary

    min_num_bad_boundary = min(pg_filtered_rt2.values())

    pg_filtered_rt3 = []

    for rt in pg_filtered_rt2.keys():
        if pg_filtered_rt2[rt] <= min_num_bad_boundary + 1:
            pg_filtered_rt3.append(rt)

    return pg_filtered_rt3


def check_if_the_top1_fragment_is_good_in_ref_sample(pg, ref_pg):

    # check the top1 fragment is a good one, otherwise, remove it

    top1 = find_top_n_fragment(1, ref_pg)

    # ref_pg_rt_apex = ref_pg['peak_rt']
    ref_pg_rt_left = ref_pg['rt_left']
    ref_pg_rt_right = ref_pg['rt_right']
    top1_rt_list = ref_pg['ms2']['rt_list'][top1]
    top1_i_list = ref_pg['ms2']['i_list'][top1]

    top1_rt_list_in_range, top1_i_list_in_range = \
        chrom.get_chrom_range(ref_pg_rt_left, ref_pg_rt_right, top1_rt_list, top1_i_list)

    if_ref_good = swath_quant.check_if_displayed_peak_a_good_one(top1_rt_list, top1_i_list, 1, 3)

    return if_ref_good


def check_if_top1_fragment_good_in_this_pg(pg, ref_pg, if_top1_fragment_good_in_ref_sample, rt):

    if_top1_fragment_good_in_this_pg = 1

    if if_top1_fragment_good_in_ref_sample == 0:
        if_top1_fragment_good_in_this_pg = 0

    else:
        top1 = find_top_n_fragment(1, ref_pg)
        top1_rt_list = pg[rt]['ms2']['rt_list'][top1]
        top1_i_list = pg[rt]['ms2']['i_list'][top1]
        ref_pg_width = ref_pg['rt_right'] - ref_pg['rt_left']
        pg_rt_left, pg_rt_right = get_boundary_for_pg_best(pg, rt, ref_pg_width)
        top1_rt_list_in_range, top1_i_list_in_range = \
            chrom.get_chrom_range(pg_rt_left, pg_rt_right, top1_rt_list, top1_i_list)
        if_top1_fragment_good_in_this_pg = swath_quant.check_if_displayed_peak_a_good_one(
            top1_rt_list, top1_i_list, 1, 3)

    return if_top1_fragment_good_in_this_pg


def compute_corr_using_fragments(pg, ref_pg, pg_filtered_rt, fragments_list, if_top1_fragment_good_in_ref_sample):

    pg_corr = {}

    for rt in pg_filtered_rt:

        if_top1_fragment_good_in_this_pg = check_if_top1_fragment_good_in_this_pg(
            pg, ref_pg, if_top1_fragment_good_in_ref_sample, rt)

        if if_top1_fragment_good_in_ref_sample == 1 and if_top1_fragment_good_in_this_pg == 0:
            fragments_list = remove_top1_fragment(ref_pg)

        x = []
        y = []

        for fragment in fragments_list:
            x.append(ref_pg['ms2']['peak_apex_i'][fragment])
            y0 = get_intensity_for_closest_rt(
                rt, pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])
            y.append(y0)

        if len(x) > 2 and len(y) == len(x):
            pg_corr[rt] = np.corrcoef(x, y)[0][1]  # note: this is R, not R2

            if math.isnan(pg_corr[rt]) == 1:
                pg_corr[rt] = -1
        else:
            pg_corr[rt] = -1

    return pg_corr


def get_list_of_rt_with_top_corr(pg_corr):

    rt_list = []
    top_corr = max(pg_corr.values())

    for rt in pg_corr.keys():
        if pg_corr[rt] >= top_corr - 0.3:
            rt_list.append(rt)

    return rt_list


def most_correlated_peak_group_based_on_fragment_intensity(pg, ref_pg, pg_filtered_rt):

    pg_filtered_rt = filter_pg_based_on_peak_boundary(pg_filtered_rt, pg)

    if len(pg_filtered_rt) > 0:

        if_top1_fragment_good_in_ref_sample = check_if_the_top1_fragment_is_good_in_ref_sample(
            pg, ref_pg)

        # select fragment list
        if if_top1_fragment_good_in_ref_sample == 1:
            fragments_list = ref_pg['ms2']['i_list'].keys()

        else:
            fragments_list = remove_top1_fragment(ref_pg)

        pg_corr = compute_corr_using_fragments(pg, ref_pg, pg_filtered_rt, fragments_list,
                                               if_top1_fragment_good_in_ref_sample)

        list_of_rt_with_top_corr = get_list_of_rt_with_top_corr(pg_corr)

        if len(list_of_rt_with_top_corr) == 1:
            return list_of_rt_with_top_corr[0]

        else:
            # if there are multiple pg with similar R2 value (in the same peak),
            # select the one with highest intensity (apex)
            top1_fragment = find_top_n_fragment(1, ref_pg)
            top2_fragment = find_top_n_fragment(2, ref_pg)
            top3_fragment = find_top_n_fragment(3, ref_pg)

            rt0 = find_best_rt_based_on_top3_fragment(
                list_of_rt_with_top_corr, pg, top1_fragment, top2_fragment, top3_fragment)

            return rt0

    else:
        print 'WARNING: no peak group found when computing fragment intensity correlation'
        return ref_pg['peak_rt']


def combine_pg_sorted_max_rt(pg_sorted_max, pg_sorted_max_spearman):

    for rt in pg_sorted_max_spearman:
        if not rt in pg_sorted_max:
            pg_sorted_max.append(rt)

    return pg_sorted_max


def remove_top1_fragment(ref_pg):

    top1 = find_top_n_fragment(1, ref_pg)

    top_list = []

    for fragment in ref_pg['ms2']['i_list'].keys():
        if fragment != top1:
            top_list.append(fragment)

    return top_list


def find_best_rt_based_on_top3_fragment(pg_sorted_max, pg, top1_fragment, top2_fragment, top3_fragment):

    # find the pg with highest intensity sum from top3 fragment
    pg_vote = {}
    for rt in pg_sorted_max:
        pg_vote[rt] = 0

    # pair wise comparison for voting
    for rt1 in pg_sorted_max:
        for rt2 in pg_sorted_max:
            if rt2 != rt1:
                vote1, vote2 = pair_wise_vote(
                    rt1, rt2, pg, top1_fragment, top2_fragment, top3_fragment)
                pg_vote[rt1] += vote1
                pg_vote[rt2] += vote2

    max_vote = max(pg_vote.values())

    pg_vote2 = {}
    for rt in pg_vote.keys():
        if pg_vote[rt] == max_vote:
            pg_vote2[rt] = pg_vote[rt]

    rt_final = -1

    if len(pg_vote2.keys()) == 1:
        rt_final = pg_vote2.keys()[0]
    else:
        rt_final = get_max_intensity_sum(pg_vote2, pg, top1_fragment, top2_fragment, top3_fragment)

    return rt_final


def get_max_intensity_sum(pg_vote2, pg, top1_fragment, top2_fragment, top3_fragment):

    rt_int_sum = {}
    for rt in pg_vote2.keys():
        i1, i2, i3 = get_intensity_for_top3_fragments(
            top1_fragment, top2_fragment, top3_fragment, rt, pg)
        rt_int_sum[rt] = i1 + i2 + i3

    max_i = max(rt_int_sum.values())

    rt_final = -1

    for rt in rt_int_sum.keys():

        if rt_int_sum[rt] == max_i:
            rt_final = rt

    return rt_final


def pair_wise_vote(rt1, rt2, pg, top1_fragment, top2_fragment, top3_fragment):

    vote1 = 0
    vote2 = 0

    i1, i2, i3 = get_intensity_for_top3_fragments(top1_fragment, top2_fragment, top3_fragment, rt1, pg)
    i1_b, i2_b, i3_b = get_intensity_for_top3_fragments(top1_fragment, top2_fragment, top3_fragment, rt2, pg)

    if i1 > i1_b:
        vote1 += 1
    else:
        vote2 += 1
    if i2 > i2_b:
        vote1 += 1
    else:
        vote2 += 1
    if i3 > i3_b:
        vote1 += 1
    else:
        vote2 += 1

    return vote1, vote2


def get_top_voted_rt(voted_rt):
    voted_rt2 = {}
    vote_max = max(voted_rt.values())
    for rt in voted_rt:
        if voted_rt[rt] == vote_max:
            voted_rt2[rt] = voted_rt[rt]
    return voted_rt2


def get_vote_based_on_top3_fragments(rt_list, top1_fragment, top2_fragment, top3_fragment, pg, i1, i2, i3):

    vote_based_on_top3_fragments = {}

    for rt in rt_list:
        i1_b, i2_b, i3_b = get_intensity_for_top3_fragments(
            top1_fragment, top2_fragment, top3_fragment, rt, pg)
        vote_based_on_top3_fragments[rt] = vote_which_rt_is_better_based_on_top3_fragments(
            i1, i2, i3, i1_b, i2_b, i3_b)

    return vote_based_on_top3_fragments


def vote_which_rt_is_better_based_on_top3_fragments(i1, i2, i3, i1_b, i2_b, i3_b):

    vote_num = 0

    if i1_b > i1:
        vote_num += 1
    if i2_b > i2:
        vote_num += 1
    if i3_b > i3:
        vote_num += 1

    return vote_num


def get_intensity_for_top3_fragments(top1_fragment, top2_fragment, top3_fragment, rt0, pg):

    i1 = get_intensity_for_closest_rt(
        rt0, pg[rt0]['ms2']['rt_list'][top1_fragment], pg[rt0]['ms2']['i_list'][top1_fragment])
    i2 = get_intensity_for_closest_rt(rt0, pg[rt0]['ms2']['rt_list'][top2_fragment],
                                      pg[rt0]['ms2']['i_list'][top2_fragment])
    i3 = get_intensity_for_closest_rt(rt0, pg[rt0]['ms2']['rt_list'][top3_fragment],
                                      pg[rt0]['ms2']['i_list'][top3_fragment])

    return i1, i2, i3


def get_max_rt_from_pg_sorted(rt_list, pg_corr):

    rt_final = []

    rt_max = rt_list[0]
    rt_final.append(rt_max)

    corr_max = pg_corr[rt_max]

    for rt in rt_list[1:]:
        if pg_corr[rt] >= corr_max - 0.3:  # empirical
            rt_final.append(rt)

    return rt_final


def get_intensity_for_closest_rt(rt0, rt_list, i_list):

    # based on a rt value, find the best match rt value, then find out the chrom peak intensity

    i = i_list[0]
    rt = rt_list[0]

    if type(rt) == float and type(rt) == float:
        rt_dif = abs(rt - rt0)

        for rt_, i_ in zip(rt_list, i_list):
            rt_dif_ = abs(rt_ - rt0)

            if rt_dif_ < rt_dif:
                i = i_
                rt = rt_
                rt_dif = rt_dif_
        return i + 0.1  # sometimes i is 0
    else:
        return 0.1


def filter_peak_group_top_fragment_peak_boundary(n, pg, ref_pg, pg_filtered_rt):

    pg_filtered_rt2 = []

    top_fragment = find_top_n_fragment(n, ref_pg)

    if top_fragment == 'NA':
        return pg_filtered_rt

    ref_sample_peak_width = float(ref_pg['rt_right'] - ref_pg['rt_left'])

    for rt in pg_filtered_rt:

        if_top_fragment_good = check_peak_group_top_fragment_peak_boundary(
            n, rt, top_fragment, pg, ref_sample_peak_width)

        if if_top_fragment_good == 1:

            pg_filtered_rt2.append(rt)  # good peak boundary

        else:
            # if this is not good, sometimes it is still the good pg. This can happen
            # due to interfering signals
            next_fragment = find_top_n_fragment(n + 1, ref_pg)

            if next_fragment != 'NA':

                if_next_fragment_good = check_peak_group_top_fragment_peak_boundary(
                    n + 1, rt, next_fragment, pg, ref_sample_peak_width)

                if if_next_fragment_good == 1:
                    pg_filtered_rt2.append(rt)

    return pg_filtered_rt2


def check_peak_group_top_fragment_peak_boundary(n, rt, fragment, pg, ref_sample_peak_width):

    if_good = 0

    # if the peak intensity is too low (below 100, it will not be a good signal in tripleTOF 5600

    if 'peak_apex_i' in pg[rt].keys():

        if fragment in pg[rt]['peak_apex_i'].keys():

            if pg[rt]['ms2']['peak_apex_i'][fragment] > 200.0:

                peak_width = float(pg[rt]['ms2']['rt_right'][fragment] - pg[rt]['ms2']['rt_left'][fragment])

                if 1.0 / parameters.PEAK_WIDTH_FOLD_VARIATION <= round(peak_width / ref_sample_peak_width, 1) <= parameters.PEAK_WIDTH_FOLD_VARIATION:
                    if_good = 1

    return if_good


def find_top_fragment_with_peak(pg):

    top_fragment = {}

    for rt in pg.keys():
        if_peak_found = 0
        top_n_fragment_found = -1
        n = 1
        for fragment in sorted(pg[rt]['ms2']['i'].keys(), reverse=True):
            for rt0 in pg[rt]['ms2']['peaks_rt'][fragment]:
                if abs(rt0 - rt) < parameters.MAX_RT_TOLERANCE:
                    if_peak_found = 1
                    break
            if if_peak_found == 1:
                top_n_fragment_found = n
                break
            else:
                n += 1
        top_fragment[rt] = top_n_fragment_found

    return top_fragment


def find_best_match_pg(pg, ref_pg, sample):

    # for debugging
    if sample == 'gold4':
        pass

    if len(pg) == 0:
        return 0

    elif len(pg) == 1:
        pg_best = only_one_pg(pg, ref_pg)
        return pg_best

    else:
        pg_best = find_best_match_pg_rule_a(pg, ref_pg, sample)
    return pg_best


def only_one_pg(pg, ref_pg):
    return get_peak_group_values(pg, pg.keys()[0], ref_pg)


def find_best_match_pg_rule_a(pg, ref_pg, sample):
    # for debugging
    if sample == 'gold40':
        pass

    # filter out peak groups without top 1 fragment as a peak
    pg_filtered_rt = filter_peak_group_top_fragment(1, pg, ref_pg, pg.keys())

    if len(pg_filtered_rt) == 1:

        pg_best = get_peak_group_values(pg, pg_filtered_rt[0], ref_pg)

    elif len(pg_filtered_rt) == 0:
        pg_best = find_best_match_pg_rule_b(pg, ref_pg, pg.keys(), sample)

    elif len(pg_filtered_rt) > 1:
        pg_best = find_best_match_pg_rule_b(pg, ref_pg, pg_filtered_rt, sample)

    return pg_best


def find_best_match_pg_rule_b(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'gold4':
        pass
    # filter out peak groups without top 2 fragment as a peak
    pg_filtered_rt2 = filter_peak_group_top_fragment(2, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_c(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_c(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_c(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'gold40':
        pass

    # filter out peak groups without MS1 as a peak
    # pg_filtered_rt2 = filter_peak_group_ms1(pg, pg_filtered_rt)
    # better not to use this because sometimes MS1 has inteference. like
    # golden set gold 40, for peptide # 531.
    pg_filtered_rt2 = pg_filtered_rt

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_d(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_d(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_f(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'BR63':
        pass

    # filter out peak groups without top 1 fragment showing good peak boundary
    pg_filtered_rt2 = filter_peak_group_top_fragment_peak_boundary(2, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_e(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_e(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_g(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'gold10':
        pass

    # filter out peak groups without top 2 fragment showing good peak boundary
    pg_filtered_rt2 = filter_peak_group_top_fragment_peak_boundary(1, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_f(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_f(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_d(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'BR1':
        pass

    # filter out peak groups without top 1 fragment showing good peak shape
    pg_filtered_rt2 = filter_peak_group_peak_shape(1, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_g(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_g(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_e(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'gold10':
        pass

    # filter out peak groups without top 2 fragment showing good peak shape
    pg_filtered_rt2 = filter_peak_group_peak_shape(2, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        pg_best = find_best_match_pg_rule_h(pg, ref_pg, pg_filtered_rt, sample)

    elif len(pg_filtered_rt2) > 1:
        pg_best = find_best_match_pg_rule_h(pg, ref_pg, pg_filtered_rt2, sample)

    return pg_best


def find_best_match_pg_rule_h(pg, ref_pg, pg_filtered_rt, sample):

    # for debugging
    if sample == 'gold11':
        pass
    # print sample  ##hahaha

    # select the peak group with highest correlation to the reference peak
    # group in terms of fragment intensity
    pg_filtered_rt2 = filter_peak_group_peak_shape(2, pg, ref_pg, pg_filtered_rt)

    if len(pg_filtered_rt2) == 1:
        pg_best = get_peak_group_values(pg, pg_filtered_rt2[0], ref_pg)

    elif len(pg_filtered_rt2) == 0:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(
            pg, ref_pg, pg_filtered_rt)
        pg_best = get_peak_group_values(pg, best_pg_rt, ref_pg)

    elif len(pg_filtered_rt2) > 1:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(
            pg, ref_pg, pg_filtered_rt2)
        pg_best = get_peak_group_values(pg, best_pg_rt, ref_pg)

    return pg_best


def build_other_sample_peak_group(chrom_data, tg, ref_pg, peak_group_candidates, sample):

    pg = data_holder.NestedDict()

    for peak_rt in peak_group_candidates[tg][sample].keys():
        # if_bad_pg = 0

        # sometimes ms1 has no signal

        if hasattr(chrom_data[tg][sample][tg], 'rt_list') \
                and hasattr(chrom_data[tg][sample][tg], 'i_list') \
                and hasattr(chrom_data[tg][sample][tg], 'peak_apex_i_list') \
                and hasattr(chrom_data[tg][sample][tg], 'peak_apex_rt_list'):

            pg[peak_rt]['ms1']['rt_list'] = chrom_data[tg][sample][tg].rt_list
            pg[peak_rt]['ms1']['i_list'] = chrom_data[tg][sample][tg].i_list
            pg[peak_rt]['ms1']['peak_apex_i_list'] = chrom_data[tg][sample][tg].peak_apex_i_list
            pg[peak_rt]['ms1']['peak_apex_rt_list'] = chrom_data[tg][sample][tg].peak_apex_rt_list

            if pg[peak_rt]['ms1']['peak_apex_rt_list'] != "NA" and pg[peak_rt]['ms1']['peak_apex_rt_list'] != "NA":
                pg[peak_rt]['ms1']['peak_apex_i'] = \
                    find_peak_apex_i_from_list(chrom_data[tg][sample][tg].peak_apex_i_list, \
                                                chrom_data[tg][sample][tg].peak_apex_rt_list, peak_rt)
            else:
                pg[peak_rt]['ms1']['peak_apex_i'] = "NA"

        else:
            pg[peak_rt]['ms1']['rt_list'] = "NA"
            pg[peak_rt]['ms1']['i_list'] = "NA"
            pg[peak_rt]['ms1']['peak_apex_i_list'] = "NA"
            pg[peak_rt]['ms1']['peak_apex_rt_list'] = "NA"
            pg[peak_rt]['ms1']['peak_apex_i'] = "NA"


        for fragment in ref_pg['ms2']['rt_list'].keys():
            if fragment != tg:
                if hasattr(chrom_data[tg][sample][fragment], 'rt_list') == 1:
                    if not chrom_data[tg][sample][fragment].rt_list == "NA":
                        pg[peak_rt]['ms2']['rt_list'][fragment] = chrom_data[
                            tg][sample][fragment].rt_list
                        pg[peak_rt]['ms2']['i_list'][fragment] = chrom_data[tg][sample][fragment].i_list
                        pg[peak_rt]['ms2']['peak_apex_i_list'][fragment] = chrom_data[
                            tg][sample][fragment].peak_apex_i_list
                        pg[peak_rt]['ms2']['peak_apex_rt_list'][fragment] = chrom_data[
                            tg][sample][fragment].peak_apex_rt_list
                        pg[peak_rt]['ms2']['peak_apex_i'][fragment] = \
                            find_peak_apex_i_from_list(chrom_data[tg][sample][fragment].peak_apex_i_list,
                                                       chrom_data[tg][sample][fragment].peak_apex_rt_list, peak_rt)

                        pg[peak_rt]['ms2']['rt_left'][fragment], pg[peak_rt]['ms2']['rt_right'][fragment] = \
                            get_peak_boundary_from_fragment(
                                peak_group_candidates, tg, sample, peak_rt, fragment)

    if len(pg) == 0:
        print 'no peak group. likely NA values'

    return pg


def get_peak_boundary_from_fragment(peak_group_candidates, tg, sample, peak_rt, fragment0):

    rt_left0 = -1.0
    rt_right0 = -1.0

    if len(peak_group_candidates[tg][sample][peak_rt].matched_fragments) > 1:

        fragments = peak_group_candidates[tg][sample][peak_rt].matched_fragments
        fragments_rt_left = peak_group_candidates[tg][
            sample][peak_rt].matched_fragments_peak_rt_left
        fragments_rt_right = peak_group_candidates[tg][
            sample][peak_rt].matched_fragments_peak_rt_right

        for fragment, rt_left, rt_right in zip(fragments, fragments_rt_left, fragments_rt_right):

            if fragment == fragment0:

                rt_left0 = rt_left
                rt_right0 = rt_right
                break

    return rt_left0, rt_right0


def build_reference_peak_group_based_on_replicates(display_data, this_unique_sample_list, chrom_data, tg):

    # select the best pg from the replicates####
    # TODO
    # NOT TO DO NOW

    ref_pg = data_holder.NestedDict()

    ref_sample = ref_sample_data[tg].sample_name
    ref_pg['peak_rt'] = ref_sample_data[tg].peak_rt_found
    ref_pg['rt_left'] = ref_sample_data[tg].peak_rt_left
    ref_pg['rt_right'] = ref_sample_data[tg].peak_rt_right

    ref_pg['ms1']['rt_list'] = chrom_data[tg][ref_sample][tg].rt_list
    ref_pg['ms1']['i_list'] = chrom_data[tg][ref_sample][tg].i_list
    ref_pg['ms1']['peak_apex_i_list'] = chrom_data[tg][ref_sample][tg].peak_apex_i_list
    ref_pg['ms1']['peak_apex_rt_list'] = chrom_data[tg][ref_sample][tg].peak_apex_rt_list
    ref_pg['ms1']['peak_apex_i'] = \
        find_peak_apex_i_from_list(chrom_data[tg][ref_sample][tg].peak_apex_i_list,
                                   chrom_data[tg][ref_sample][tg].peak_apex_rt_list, ref_pg['peak_rt'])

    display_data[tg][ref_sample]['ms1']['peak_apex_i'] = ref_pg['ms1']['peak_apex_i']

    for fragment in chrom_data[tg][ref_sample].keys():
        if fragment != tg:
            ref_pg['ms2']['rt_list'][fragment] = chrom_data[tg][ref_sample][fragment].rt_list
            ref_pg['ms2']['i_list'][fragment] = chrom_data[tg][ref_sample][fragment].i_list
            ref_pg['ms2']['peak_apex_i_list'][fragment] = chrom_data[
                tg][ref_sample][fragment].peak_apex_i_list
            ref_pg['ms2']['peak_apex_rt_list'][fragment] = chrom_data[
                tg][ref_sample][fragment].peak_apex_rt_list
            ref_pg['ms2']['peak_apex_i'][fragment] = find_peak_apex_i_from_list(chrom_data[tg][ref_sample][fragment].peak_apex_i_list,
                                                                                chrom_data[tg][ref_sample][fragment].peak_apex_rt_list, ref_pg['peak_rt'])

            display_data[tg][ref_sample]['ms2']['peak_apex_i'][
                fragment] = ref_pg['ms2']['peak_apex_i'][fragment]

    return display_data, ref_pg


def build_reference_peak_group(display_data, ref_sample_data, chrom_data, tg):

    ref_pg = data_holder.NestedDict()

    ref_sample = ref_sample_data[tg].sample_name
    ref_pg['peak_rt'] = ref_sample_data[tg].peak_rt_found
    ref_pg['rt_left'] = ref_sample_data[tg].peak_rt_left
    ref_pg['rt_right'] = ref_sample_data[tg].peak_rt_right

    # sometimes MS1 of the reference sample is NA

    if hasattr(chrom_data[tg][ref_sample][tg], 'rt_list') \
        and hasattr(chrom_data[tg][ref_sample][tg], 'i_list') \
        and hasattr(chrom_data[tg][ref_sample][tg], 'peak_apex_i_list') \
        and hasattr(chrom_data[tg][ref_sample][tg], 'peak_apex_rt_list'):

        ref_pg['ms1']['rt_list'] = chrom_data[tg][ref_sample][tg].rt_list
        ref_pg['ms1']['i_list'] = chrom_data[tg][ref_sample][tg].i_list
        ref_pg['ms1']['peak_apex_i_list'] = chrom_data[tg][ref_sample][tg].peak_apex_i_list
        ref_pg['ms1']['peak_apex_rt_list'] = chrom_data[tg][ref_sample][tg].peak_apex_rt_list

        if chrom_data[tg][ref_sample][tg].rt_list == "NA":
            ref_pg['ms1']['peak_apex_i'] = "NA"
        else:
            ref_pg['ms1']['peak_apex_i'] = \
                find_peak_apex_i_from_list(chrom_data[tg][ref_sample][tg].peak_apex_i_list,
                                               chrom_data[tg][ref_sample][tg].peak_apex_rt_list, ref_pg['peak_rt'])

        display_data[tg][ref_sample]['ms1']['peak_apex_i'] = ref_pg['ms1']['peak_apex_i']

        for fragment in chrom_data[tg][ref_sample].keys():
            if fragment != tg:
                ref_pg['ms2']['rt_list'][fragment] = chrom_data[tg][ref_sample][fragment].rt_list
                ref_pg['ms2']['i_list'][fragment] = chrom_data[tg][ref_sample][fragment].i_list
                ref_pg['ms2']['peak_apex_i_list'][fragment] = chrom_data[
                    tg][ref_sample][fragment].peak_apex_i_list
                ref_pg['ms2']['peak_apex_rt_list'][fragment] = chrom_data[
                    tg][ref_sample][fragment].peak_apex_rt_list
                ref_pg['ms2']['peak_apex_i'][fragment] = find_peak_apex_i_from_list(chrom_data[tg][ref_sample][fragment].peak_apex_i_list,
                                                                                    chrom_data[tg][ref_sample][fragment].peak_apex_rt_list, ref_pg['peak_rt'])

                display_data[tg][ref_sample]['ms2']['peak_apex_i'][
                    fragment] = ref_pg['ms2']['peak_apex_i'][fragment]

    else:
        ref_pg = "NA"

    return display_data, ref_pg

def find_peak_apex_i_from_list(peak_apex_i, peak_apex_rt, peak_rt):
    if len(peak_apex_i) == 1:
        return peak_apex_i[0]

    else:
        rt0 = peak_apex_rt[0]
        i0 = peak_apex_i[0]
        rt_dif = abs(peak_rt - rt0)
        for rt, i in zip(peak_apex_rt, peak_apex_i):
            rt_dif2 = abs(rt - peak_rt)
            if rt_dif2 < rt_dif:
                rt0 = rt
                i0 = i
                rt_dif = rt_dif2

        return i0


def find_all_rt_values(chrom_data, tg, sample):

    all_rt = []
    for fragment in chrom_data[tg][sample].keys():
        if hasattr(chrom_data[tg][sample][fragment], 'peak_apex_rt_list'):
            peak_apex_rt_list = chrom_data[tg][sample][fragment].peak_apex_rt_list
            if peak_apex_rt_list != 'NA':
                all_rt += peak_apex_rt_list  # include every element in the 2nd list

    all_rt = sorted(list(set(all_rt)))

    for i in range(10):  # repeat the following 10 times
        all_rt = binning_rt_values(all_rt)

    return all_rt


def find_peak_group_candidates(chrom_data, sample_id):

    peak_group_candidates = data_holder.NestedDict()

    for tg in chrom_data.keys():

        for sample in sample_id:

            # if sample == 'gold10':
            #     pass

            all_rt = find_all_rt_values(chrom_data, tg, sample)

            if len(all_rt) == 0:
                # no peak found, this sample has no rt. Most likely this sample has no chrom.
                continue
            else:

                for rt in all_rt:

                    # compute the peak boundary for each fragment, not the consensus peak boundary

                    # if sample == 'gold10' and abs(rt - 3893.9) < 1:
                    #     pass

                    this_peak_group = data_holder.PeakGroup(chrom_data, tg, sample, rt)

                    if this_peak_group.num_matched_fragments >= 3:
                        peak_group_candidates[tg][sample][rt] = this_peak_group

                if len(peak_group_candidates[tg][sample].keys()) == 0:
                    # if no good peak group is found, then use whatever peak group
                    for rt in all_rt:

                        # compute the peak boundary for each fragment, not the consensus peak
                        # boundary
                        this_peak_group = data_holder.PeakGroup(chrom_data, tg, sample, rt)
                        peak_group_candidates[tg][sample][rt] = this_peak_group

                    if len(peak_group_candidates[tg][sample].keys()) == 0:
                        # if still no peak group is found. Most likely in case of empty chrom. Use
                        # all rt as peak group
                        for rt in all_rt:
                            this_peak_group = data_holder.PeakGroup(chrom_data, tg, sample, rt)
                            peak_group_candidates[tg][sample][rt] = this_peak_group

    return peak_group_candidates


def get_max_num_good_fragments(peak_group_candidates, tg, sample):

    n = []
    for rt in peak_group_candidates[tg][sample].keys():
        n1 = len(peak_group_candidates[tg][sample][rt].matched_fragments)
        n.append(n1)
    if len(n) > 0:
        return max(n)
    else:
        return 0


def check_best_peak_group_from_reference_sample(ref_sample_data, peak_group_candidates, tg, chrom_data):

    # find out the peak group with highest number of good fragments

    sample = ref_sample_data[tg].sample_name

    # go through all rt for pg in the peak_group_candidates,
    # find the peak groups with highest number of fragments.
    # also get the MS1 status
    good_fragments = []
    if_ms1 = 0
    rt_dif = parameters.MAX_RT_TOLERANCE
    rt_found = -1.0

    selected_rt_num_fragment = {}
    selected_rt_if_ms1 = {}

    max_num_good_fragments = get_max_num_good_fragments(peak_group_candidates, tg, sample)

    for rt in peak_group_candidates[tg][sample].keys():

        num_fragments = len(peak_group_candidates[tg][sample][rt].matched_fragments)

        if num_fragments >= max_num_good_fragments - 3:

            good_fragments = peak_group_candidates[tg][sample][rt].matched_fragments
            rt_dif = abs(rt - ref_sample_data[tg].peak_rt)
            rt_found = rt

            selected_rt_num_fragment[rt] = num_fragments
            selected_rt_if_ms1[rt] = peak_group_candidates[tg][sample][rt].if_ms1_peak

    # sometimes there are multiple rt with max number of good fragment
    # in this case, check ms1, peak shape
    # case: golden data set #96. gold80

    selected_rt_num_fragment_2, selected_rt_if_ms1_2 = \
        select_peak_with_high_number_of_good_fragment(
            selected_rt_num_fragment, selected_rt_if_ms1, max_num_good_fragments)

    if len(selected_rt_num_fragment_2.keys()) > 1:

        selected_rt_num_fragment_3, selected_rt_if_ms1_3 = select_rt_from_reference_sample_based_on_ms1(
            selected_rt_num_fragment_2, selected_rt_if_ms1_2)

        if len(selected_rt_num_fragment_3.keys()) == 1:
            rt_found = selected_rt_num_fragment_3.keys()[0]
            num_good_fragments = selected_rt_num_fragment_3[rt_found]
            if_ms1 = selected_rt_if_ms1_3[rt_found]
            good_fragments = peak_group_candidates[tg][sample][rt_found].matched_fragments
            rt_dif = abs(rt_found - ref_sample_data[tg].peak_rt)
            return float(rt_found), num_good_fragments, if_ms1, good_fragments, float(rt_dif)

        elif len(selected_rt_num_fragment_3.keys()) > 1:
            # select based on ms1 and ms2 peak shape
            num_good_shape_peak_ms2_dict = {}
            for rt in selected_rt_num_fragment_3.keys():
                num_good_shape_peak_ms2_dict[rt] = get_good_shape_peak_ms2(
                    rt, sample, tg, chrom_data, peak_group_candidates)

            selected_rt_num_fragment_4 = rank_good_shape_peaks_ms2(num_good_shape_peak_ms2_dict)

            # this may be further improved in the future
            # consider intensity ranking, et al.

            rt_found = selected_rt_num_fragment_4.keys()[0]
            num_good_fragments = selected_rt_num_fragment_4[rt_found]
            if_ms1 = selected_rt_if_ms1_3[rt_found]
            good_fragments = peak_group_candidates[tg][sample][rt_found].matched_fragments
            rt_dif = abs(rt_found - ref_sample_data[tg].peak_rt)
            return float(rt_found), num_good_fragments, if_ms1, good_fragments, float(rt_dif)

    elif len(selected_rt_num_fragment_2.keys()) == 1:
        # if only one
        return float(rt_found), max_num_good_fragments, selected_rt_if_ms1_2.values()[0], good_fragments, float(rt_dif)
    else:
        # can not happen
        print "error: can not find a good rt, this should not happen. Remove this tg ", tg
        return "NA", "NA", "NA", "NA", "NA"


def rank_good_shape_peaks_ms2(good_shape_fragment):
    good_shape_fragment2 = {}
    max_num_fragments = max(good_shape_fragment.values())
    for rt in good_shape_fragment.keys():
        if good_shape_fragment[rt] == max_num_fragments:
            good_shape_fragment2[rt] = good_shape_fragment[rt]
    return good_shape_fragment2


def get_good_shape_peak_ms2(rt, sample, tg, chrom_data, peak_group_candidates):

    good_shape_fragments = []
    # get the peak boundary for the reference sample
    fragments = peak_group_candidates[tg][sample][rt].matched_fragments
    i = peak_group_candidates[tg][sample][rt].matched_fragments_i
    rt_left_list = peak_group_candidates[tg][sample][rt].matched_fragments_peak_rt_left
    rt_right_list = peak_group_candidates[tg][sample][rt].matched_fragments_peak_rt_right

    if len(fragments) > 0 and len(i) >0 and len(rt_left_list) > 0  and len(rt_right_list) > 0:

        ref_sample_rt_left, ref_sample_rt_right = chrom.get_peak_group_boundary(
            fragments, i, rt_left_list, rt_right_list, rt)

        for fragment in fragments:

            peak_intensity_apex = get_intensity_for_closest_rt(
                rt, chrom_data[tg][sample][fragment].rt_list, chrom_data[tg][sample][fragment].i_list)
            peak_intensity_left = get_intensity_for_closest_rt(
                ref_sample_rt_left, chrom_data[tg][sample][fragment].rt_list, chrom_data[tg][sample][fragment].i_list)
            peak_intensity_right = get_intensity_for_closest_rt(
                ref_sample_rt_right, chrom_data[tg][sample][fragment].rt_list, chrom_data[tg][sample][fragment].i_list)

            fold_change_left = float(peak_intensity_apex) / (float(peak_intensity_left) + 1)
            fold_change_right = float(peak_intensity_apex) / (float(peak_intensity_right) + 1)

            if fold_change_left > parameters.PEAK_SHAPE_FOLD_VARIATION and fold_change_right > parameters.PEAK_SHAPE_FOLD_VARIATION:

                good_shape_fragments.append(fragment)

    return len(good_shape_fragments)


def select_rt_from_reference_sample_based_on_ms1(selected_rt_num_fragment, selected_rt_if_ms1):

    selected_rt_num_fragment_2 = {}
    selected_rt_if_ms1_2 = {}

    for rt in selected_rt_num_fragment.keys():
        if selected_rt_if_ms1[rt] == 1:
            selected_rt_num_fragment_2[rt] = selected_rt_num_fragment[rt]
            selected_rt_if_ms1_2[rt] = selected_rt_if_ms1[rt]

    if len(selected_rt_num_fragment_2.keys()) > 0:
        return selected_rt_num_fragment_2, selected_rt_if_ms1_2
    else:
        return selected_rt_num_fragment, selected_rt_if_ms1


def select_peak_with_high_number_of_good_fragment(selected_rt_num_fragment, selected_rt_if_ms1, max_good_fragments):

    selected_rt_num_fragment_2 = {}
    selected_rt_if_ms1_2 = {}

    for rt in selected_rt_num_fragment.keys():
        if selected_rt_num_fragment[rt] >= max_good_fragments - 3:
            selected_rt_num_fragment_2[rt] = selected_rt_num_fragment[rt]
            selected_rt_if_ms1_2[rt] = selected_rt_if_ms1[rt]

    return selected_rt_num_fragment_2, selected_rt_if_ms1_2


def find_rt_for_reference_sample(ref_sample_data, peak_group_candidates, tg, chrom_data):

    # 2015.11.18. sometimes OpenSWATH picked the wrong peak group even for the reference sample.
    # here implement a check of correct peak group

    best_rt_after_check, num_good_fragments_check, if_ms1_check, good_fragments_check, rt_dif_check = \
        check_best_peak_group_from_reference_sample(
            ref_sample_data, peak_group_candidates, tg, chrom_data)

    if best_rt_after_check == "NA":

        return "NA", "NA", "NA"

    else:

        return good_fragments_check, rt_dif_check, best_rt_after_check


def find_closest_rt_value(rt0, rt_list):
    rt1 = rt_list[0]
    rt_dif = 9999.99
    for rt in rt_list[1:]:
        rt_dif2 = abs(rt1 - rt)
        if rt_dif2 < rt_dif:
            rt1 = rt
            rt_dif = rt_dif2
    return rt1


def refine_peak_forming_fragments_based_on_reference_sample(ref_sample_data, chrom_data, peptide_data, peak_group_candidates):

    # for each tg, find the peak_forming fragments in the reference sample

    for tg in chrom_data.keys():

        print tg, "refine_peak_forming_fragments_based_on_reference_sample"

        ref_sample = ref_sample_data[tg].sample_name

        # find the matched rt in the reference sample
        good_fragments, rt_dif, rt_found = find_rt_for_reference_sample(
            ref_sample_data, peak_group_candidates, tg, chrom_data)

        if good_fragments == "NA" and rt_found == "NA":

            ref_sample_data.pop(tg, None)

            continue

        ref_sample_data[tg].read_peak_rt_found(rt_found)

        # if num of good fragments < MIN_FRAGMENTS, remove the tg from all variables
        if len(good_fragments) < parameters.MIN_FRAGMENTS:
            del chrom_data[tg]
            del peptide_data[tg]
            del ref_sample_data[tg]
            del peak_group_candidates[tg]
        else:
            # remove rt for other peak groups for the reference sample, only keep the
            # picked peak group rt
            for rt in peak_group_candidates[tg][ref_sample].keys():
                if rt != rt_found:
                    del peak_group_candidates[tg][ref_sample][rt]

            # for all samples, delete chrom data for non-selected fragments
            for sample in chrom_data[tg].keys():
                for fragment in chrom_data[tg][sample].keys():
                    if fragment in good_fragments or fragment == tg:
                        pass
                    else:
                        if fragment in chrom_data[tg][sample].keys():
                            del chrom_data[tg][sample][fragment]

                        if fragment in peptide_data[tg]['ms2'].keys():
                            del peptide_data[tg]['ms2'][fragment]

    return ref_sample_data, chrom_data, peptide_data, peak_group_candidates


def binning_rt_values(rt_list):

    if len(rt_list) <= 1:
        return rt_list

    rt_list2 = sorted(set(rt_list))
    rt_list3 = []
    if_skip = 0
    for i in range(1, len(rt_list2)):
        if if_skip == 1:
            if_skip = 0
            continue
        else:
            if (rt_list2[i] - rt_list2[i - 1]) > parameters.BINNING_RT_VALUE_TOLERANCE:
                rt_list3.append(rt_list2[i - 1])
            else:
                new_rt = 0.5 * (rt_list2[i] + rt_list2[i - 1])
                rt_list3.append(new_rt)
                if_skip = 1

    # the last element is missing, process it separately
    if rt_list2[-1] - rt_list3[-1] > parameters.BINNING_RT_VALUE_TOLERANCE:
        rt_list3.append(rt_list2[-1])
    else:
        new_rt = 0.5 * (rt_list2[-1] + rt_list3[-1])
        rt_list3[-1] = new_rt
    return rt_list3
