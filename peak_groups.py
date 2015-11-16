__author__ = 'guot'

import chrom
import peaks
import numpy as np
from collections import defaultdict
import parameters
import data_holder



def check_if_ms1_peak(chrom_data, tg, sample, rt):

    if_ms1_peak = 0

    peak_apex_rt = chrom_data[tg][sample][tg].peak_apex_rt
    peak_apex_i = chrom_data[tg][sample][tg].peak_apex_i
    if peak_apex_rt != 'NA':
        for rt0, i in zip(peak_apex_rt, peak_apex_i):
            if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE : # 10 s for 2hr gradient
                if_ms1_peak = 1
                break

    return if_ms1_peak


def find_matched_fragments(chrom_data, tg, sample, rt):

    matched_fragments = []
    matched_fragments_rt = []
    matched_fragments_i = []
    matched_fragments_peak_rt_left = []
    matched_fragments_peak_rt_right = []

    for fragment in chrom_data[tg][sample].keys():
        if fragment != tg:
            peak_apex_rt = chrom_data[tg][sample][fragment].peak_apex_rt
            peak_apex_i = chrom_data[tg][sample][fragment].peak_apex_i
            if peak_apex_rt != 'NA':
                for rt0, i in zip(peak_apex_rt, peak_apex_i):
                    if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE : # 10 s for 2hr gradient
                        matched_fragments.append(fragment)
                        matched_fragments_rt.append(rt0)
                        matched_fragments_i.append(i)
                        rt_list = chrom_data[tg][sample][fragment].rt_list
                        i_list = chrom_data[tg][sample][fragment].i_list
                        rt_left, rt_right = chrom.get_peak_boundary(rt_list, i_list, rt0)
                        matched_fragments_peak_rt_left.append(rt_left)
                        matched_fragments_peak_rt_right.append(rt_right)

                        break
    return matched_fragments, matched_fragments_rt, matched_fragments_i, matched_fragments_peak_rt_left, matched_fragments_peak_rt_right

def find_best_peak_group_based_on_reference_sample(display_data, ref_sample_data, chrom_data, peptide_data, peak_group_candidates, sample_id):

    for tg in chrom_data.keys():
        # build a ngram object for the peak group from reference sample
        ref_pg = build_reference_peak_group(ref_sample_data, chrom_data, tg)

        for sample in sample_id:
            if sample != ref_sample_data[tg].sample_name:
                # for each peak group, create a data structure containing all the information above
                pg = build_other_sample_peak_group(chrom_data, tg, ref_pg, peak_group_candidates, sample)

                pg_best = find_best_match_pg(pg, ref_pg)

                if pg_best == 0:
                    #TODO no peak found. cannot be!!!
                    ##  debug debug debug......

                # write into display_pg
                display_data[tg][sample]['rt_left'] = pg_best['rt_left']
                display_data[tg][sample]['rt_right'] = pg_best['rt_right']

                for fragment in pg_best['ms2']['rt_list'].keys():
                    display_data[tg][sample]['ms2']['rt_list'][fragment] = pg_best['ms2']['rt_list'][fragment]
                    display_data[tg][sample]['ms2']['i_list'][fragment] = pg_best['ms2']['i_list'][fragment]
                    display_data[tg][sample]['ms2']['peak_apex_i'][fragment] = pg_best['ms2']['peak_apex_i'][fragment]
                    display_data[tg][sample]['ms2']['if_found_peak'][fragment] = pg_best['ms2']['if_found_peak'][fragment]

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

def get_peak_group_values(pg, rt, ref_pg):

    pg_best = data_holder.Nested_dict()

    pg_best['peak_rt'] = rt

    # use peak width from ref_pg to compute the boundary of this peak group
    pg_best['rt_left'] = rt - 0.5 * (ref_pg['rt_right'] - ref_pg['rt_left'])
    pg_best['rt_right'] = rt + 0.5 * (ref_pg['rt_right'] - ref_pg['rt_left'])

    for fragment in pg[rt]['ms2']['rt_list'].keys():

        # get the rt and i list only for the peak group
        pg_best['ms2']['rt_list'][fragment], pg_best['ms2']['i_list'][fragment] = \
            chrom.get_chrom_range(pg_best['rt_left'], pg_best['rt_right'],
                                  pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])

        pg_best['ms2']['peak_apex_i'][fragment], pg_best['ms2']['if_found_peak'][fragment] = get_fragment_intensity_for_peak_group(\
            pg[rt]['ms2']['peak_apex_rt'][fragment], pg[rt]['ms2']['peak_apex_i'][fragment], rt)

    # get the rt and i list only for the peak group
    pg_best['ms1']['rt_list'], pg_best['ms1']['i_list'] = \
        chrom.get_chrom_range(pg_best['rt_left'], pg_best['rt_right'], pg[rt]['ms1']['rt_list'], pg[rt]['ms1']['i_list'])

    pg_best['ms1']['i'], pg_best['ms1']['if_found_peak'] = get_fragment_intensity_for_peak_group(\
            pg[rt]['ms1']['peak_apex_rt'], pg[rt]['ms1']['peak_apex_i'], rt)

    return pg_best

def find_top_n_fragment(option, rt, pg):
    # sort fragment based on i, and then select the top n fragment
    fragment_sorted = []

    for fragment, i in sorted(zip(pg[rt]['ms2']['i'].values(), pg[rt]['ms2']['i'].values()), reverse = True):
        fragment_sorted.append(fragment)

    return fragment_sorted[int(option) - 1]

def filter_peak_group_top_fragment(n, pg):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        if_peak_found = 0
        for rt0 in pg[rt]['ms2']['peaks_rt'][fragment]:
            if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:
                if_peak_found = 1
                break
        if if_peak_found == 0:
            del pg[rt]
    return pg

def filter_peak_group_ms1(pg):
    for rt in pg.keys():
        if_peak_found = 0
        for rt0 in pg[rt]['ms1']['peaks_rt']:
            if abs(rt - rt0) < parameters.MAX_RT_TOLERANCE:
                if_peak_found = 1
                break
        if if_peak_found == 0:
            del pg[rt]
    return pg

def compute_transition_intensity(rt0, rt_list, i_list):
    #based on a rt value, find the best match rt value, then find out the chrom peak intensity
    transition_intensity = {}
    for this_transition in rt_list.keys():
        rt_distance = 9999
        this_intensity = 0.0
        for rt, i in zip (rt_list[this_transition], i_list[this_transition]):
            if abs(rt - rt0) < rt_distance:
                rt_distance = abs(rt - rt0)
                this_intensity = i
        transition_intensity[this_transition] = this_intensity
    return transition_intensity

def filter_peak_group_peak_shape(n, pg):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        peak_intensity_apex = compute_transition_intensity(rt, pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        peak_intensity_left = compute_transition_intensity(pg[rt]['ms2']['rt_left'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        peak_intensity_right = compute_transition_intensity(pg[rt]['ms2']['rt_right'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        fold_change_left = float(peak_intensity_apex) / float(peak_intensity_left)
        fold_change_right = float(peak_intensity_apex) / float(peak_intensity_right)

        if  fold_change_left > parameters.PEAK_SHAPE_FOLD_VARIATION and fold_change_right > parameters.PEAK_SHAPE_FOLD_VARIATION:
            pass  # good peak boundary
        else:
            del pg[rt]
    return pg

def most_correlated_peak_group_based_on_fragment_intensity(pg, ref_pg):

    pg_corr = {}

    for rt in pg.keys():

        x = [ref_pg['ms2']['peak_apex_i'][fragment] for fragment in pg[rt]['ms2']['peak_apex_i'].keys()]
        y = []
        for fragment in pg[rt]['ms2']['peak_apex_i'].keys():
            y0 = get_intensity_for_closest_rt(pg[rt]['ms2']['peak_apex_i'][fragment], pg[rt]['ms2']['peak_apex_rt'][fragment], rt)
            y.append(y0)
        pg_corr[rt] = np.corrcoef(x,y)[0][1]   #note: this is R, not R2

    # select the peak group with highest corr value
    pg_sorted = [rt for cor, rt in sorted(zip(pg_corr.values(), pg_corr.keys()), reverse=True)]

    return pg_sorted[0]

def get_intensity_for_closest_rt(peaks_i, peaks_rt, rt):
    i = -1
    rt_dif = 999.0
    for i0, rt0 in zip(peaks_i, peaks_rt):
        if abs(rt0 - rt) < rt_dif:
            rt_dif = abs(rt0 - rt)
            i = i0
    return i

def filter_peak_group_top_fragment_peak_boundary(n, pg, ref_pg):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        ref_sample_peak_width = float(ref_pg['rt_right'] - ref_pg['rt_left'])
        peak_width = float(pg[rt]['ms2']['rt_right'][fragment] - pg[rt]['ms2']['rt_left'][fragment])
        if 1.0 / parameters.PEAK_WIDTH_FOLD_VARIATION < peak_width / ref_sample_peak_width < parameters.PEAK_WIDTH_FOLD_VARIATION:
            pass  # good peak boundary
        else:
            del pg[rt]
    return pg

def find_top_fragment_with_peak(pg):

    top_fragment = {}

    for rt in pg.keys():
        if_peak_found = 0
        top_n_fragment_found = -1
        n = 1
        for fragment in sorted(pg[rt]['ms2']['i'].keys(), reverse = True):
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



def find_best_match_pg(pg, ref_pg):

    if len(pg) == 0:
        return 0
    elif len(pg) == 1:
        pg_best = only_one_pg(pg, ref_pg)
        return pg_best
    else:
        pg_best = find_best_match_pg_rule_a(pg, ref_pg)
        return pg_best

def only_one_pg(pg, ref_pg):
    return get_peak_group_values(pg, pg.keys()[0], ref_pg)


def find_best_match_pg_rule_a(pg, ref_pg):

    # filter out peak groups without top 1 fragment as a peak
    pg2 = filter_peak_group_top_fragment(1, pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_b(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_b(pg2, ref_pg)

    return pg_best


def find_best_match_pg_rule_b(pg, ref_pg):

    # filter out peak groups without top 2 fragment as a peak
    pg2 = filter_peak_group_top_fragment(2, pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_c(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_c(pg2, ref_pg)

    return pg_best

def find_best_match_pg_rule_c(pg, ref_pg):

    # filter out peak groups without MS1 as a peak
    pg2 = filter_peak_group_ms1(pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_d(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_d(pg2, ref_pg)

    return pg_best


def find_best_match_pg_rule_d(pg, ref_pg):

    # filter out peak groups without top 1 fragment showing good peak boundary
    pg2 = filter_peak_group_top_fragment_peak_boundary(1, pg, ref_pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_e(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_e(pg2, ref_pg)

    return pg_best

def find_best_match_pg_rule_e(pg, ref_pg):

    # filter out peak groups without top 2 fragment showing good peak boundary
    pg2 = filter_peak_group_top_fragment_peak_boundary(2, pg, ref_pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_f(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_f(pg2, ref_pg)

    return pg_best

def find_best_match_pg_rule_f(pg, ref_pg):

    # filter out peak groups without top 1 fragment showing good peak shape
    pg2 = filter_peak_group_peak_shape(1, pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_g(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_g(pg2, ref_pg)

    return pg_best

def find_best_match_pg_rule_g(pg, ref_pg):

    # filter out peak groups without top 2 fragment showing good peak shape
    pg2 = filter_peak_group_peak_shape(2, pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_h(pg, ref_pg)

    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_h(pg2, ref_pg)

    return pg_best

def find_best_match_pg_rule_h(pg, ref_pg):

    # select the peak group with highest correlation to the reference peak group in terms of fragment intensity
    pg2 = filter_peak_group_peak_shape(2, pg)

    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)

    elif len(pg2) == 0:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(pg, ref_pg)
        pg_best = get_peak_group_values(pg, best_pg_rt, ref_pg)

    elif len(pg2) > 1:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(pg2, ref_pg)
        pg_best = get_peak_group_values(pg2, best_pg_rt, ref_pg)

    return pg_best


def build_other_sample_peak_group(chrom_data, tg, ref_pg, peak_group_candidates, sample):

    pg = data_holder.Nested_dict()

    for peak_rt in peak_group_candidates[tg][sample].keys():
        pg[peak_rt]['ms1']['rt_list'] = chrom_data[tg][sample][tg].rt_list
        pg[peak_rt]['ms1']['i_list'] = chrom_data[tg][sample][tg].i_list
        pg[peak_rt]['ms1']['peak_apex_rt'] = chrom_data[tg][sample][tg].peak_apex_rt
        pg[peak_rt]['ms1']['peak_apex_i'] = chrom_data[tg][sample][tg].peak_apex_i

        for fragment in ref_pg['ms2']['rt_list'].keys():
            if fragment != tg:
                pg[peak_rt]['ms2']['rt_list'][fragment] = chrom_data[tg][sample][fragment]['rt_list']
                pg[peak_rt]['ms2']['i_list'][fragment] = chrom_data[tg][sample][fragment]['i_list']
                pg[peak_rt]['ms2']['peak_apex_rt'][fragment] = chrom_data[tg][sample][fragment]['peak_apex_rt']
                pg[peak_rt]['ms2']['peak_apex_i'][fragment] = chrom_data[tg][sample][fragment]['peak_apex_i']
                pg[peak_rt]['ms2']['rt_left'][fragment] = peak_group_candidates[tg][sample][peak_rt].matched_fragments_rt_left
                pg[peak_rt]['ms2']['rt_right'][fragment] = peak_group_candidates[tg][sample][peak_rt].matched_fragments_rt_right

    return pg

def build_reference_peak_group(ref_sample_data, chrom_data, tg):

    ref_pg = data_holder.Nested_dict()

    ref_sample = ref_sample_data[tg].sample_name
    ref_pg['peak_rt'] = ref_sample_data[tg].peak_rt_found
    ref_pg['rt_left'] = ref_sample_data[tg].peak_rt_left
    ref_pg['rt_right'] = ref_sample_data[tg].peak_rt_right

    ref_pg['ms1']['rt_list'] = chrom_data[tg][ref_sample][tg].rt_list
    ref_pg['ms1']['i_list'] = chrom_data[tg][ref_sample][tg].i_list
    ref_pg['ms1']['peak_apex_i'] = chrom_data[tg][ref_sample][tg].peak_apex_i

    for fragment in chrom_data[tg][ref_sample].keys():
        if fragment != tg:
            ref_pg['ms2']['rt_list'][fragment] = chrom_data[tg][ref_sample][fragment].rt_list
            ref_pg['ms2']['i_list'][fragment] = chrom_data[tg][ref_sample][fragment].i_list
            ref_pg['ms2']['peak_apex_i'][fragment] = chrom_data[tg][ref_sample][fragment].peak_apex_i

    return ref_pg

def find_all_rt_values(chrom_data, tg, sample):

    all_rt = []
    for fragment in chrom_data[tg][sample].keys():
        peak_apex_rt = chrom_data[tg][sample][fragment].peak_apex_rt
        if peak_apex_rt != 'NA':
            all_rt += peak_apex_rt # incluce every element in the 2nd list

    all_rt = sorted(list(set(all_rt)))

    for i in range(10):  #repeat the following 10 times
        all_rt = binning_rt_values(all_rt)

    return all_rt


def find_peak_group_candidates(chrom_data, sample_id):

    peak_group_candidates = data_holder.Nested_dict()

    for tg in chrom_data.keys():

        for sample in sample_id:

            all_rt = find_all_rt_values(chrom_data, tg, sample)

            if len(all_rt) == 0:
                # no peak found
                break
            else:

                for rt in all_rt:

                    #compute the peak boundary for each fragment, not the consensus peak boundary
                    this_peak_group = data_holder.Peak_group(chrom_data, tg, sample, rt)
                    if this_peak_group.num_matched_fragments >= parameters.MIN_FRAGMENTS:
                        peak_group_candidates[tg][sample][rt] = this_peak_group

    return peak_group_candidates


def find_rt_for_reference_sample(ref_sample_data, peak_group_candidates, tg):

    # maybe multiple peak groups are found for the reference sample, here, find the peak rt closest to the rt found by openswath
    sample = ref_sample_data[tg].sample_name
    best_rt = ref_sample_data[tg].peak_rt
    num_good_fragments = 0
    good_fragments = []
    rt_dif = parameters.MAX_RT_TOLERANCE
    rt_found = -1.0
    for rt in peak_group_candidates[tg][sample].keys():
        if abs(rt - best_rt) < rt_dif:
            num_fragments = len(peak_group_candidates[tg][sample][rt].matched_fragments)
            if num_fragments > num_good_fragments:
                good_fragments = peak_group_candidates[tg][sample][rt].matched_fragments
                num_good_fragments = len(good_fragments)
                rt_dif = abs(rt - ref_sample_data[tg].peak_rt)
                rt_found = rt

    return good_fragments, rt_dif, rt_found


def refine_peak_forming_fragments_based_on_reference_sample(ref_sample_data, chrom_data, peptide_data, peak_group_candidates):

    # for each tg, find the peak_forming fragments in the reference sample

    for tg in chrom_data.keys():

        ref_sample = ref_sample_data[tg].sample_name

        # find the matched rt in the reference sample
        good_fragments, rt_dif, rt_found = find_rt_for_reference_sample(ref_sample_data, peak_group_candidates, tg)

        ref_sample_data[tg].read_peak_rt_found(rt_found)

        # if num of good fragments < MIN_FRAGMENTS, remove the tg from all variables
        if len(good_fragments) < parameters.MIN_FRAGMENTS:
            del chrom_data[tg]
            del peptide_data[tg]
            del ref_sample_data[tg]
            del peak_group_candidates[tg]
        else:
            # remove rt for other peak groups for the reference sample, only keep the picked peak group rt
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

def binning_rt_values (rt):
    if len(rt) > 1:
        rt2 = sorted(list(set(rt)))
        rt3 = []
        if_skip = 0
        for i in range(1, len(rt2)):
            if if_skip == 1:
                if_skip = 0
                continue
            else:
                if (rt2[i] - rt2[i-1]) > 10 :
                    rt3.append(rt2[i-1])
                else:
                    new_rt = 0.5 * (rt2[i] + rt2[i-1])
                    rt3.append(new_rt)
                    if_skip = 1

        #the last element is missing, process it separately
        if rt2[-1] - rt3[-1] >10:
            rt3.append(rt2[-1])
        else:
            new_rt = 0.5 * (rt2[-1] + rt3[-1])
            rt3.append(new_rt)
            del rt3[-2]
        return rt3
    else:
        return rt