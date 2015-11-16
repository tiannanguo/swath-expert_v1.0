__author__ = 'guot'

import chrom
import peaks
import numpy as np
from collections import defaultdict
import parameters

import data_holder


class Peak_group(object):
    def __init__(self, all_rt, chrom_data, peptide_data, tg, sample, rt):
        self.rt = rt
        self.matched_fragments, self.matched_fragments_rt, self.matched_fragments_i = find_matched_fragments(chrom_data, peptide_data, tg, sample, rt)
        self.if_ms1_peak = check_if_ms1_peak()

def find_matched_fragments(chrom_data, peptide_data, tg, sample, rt):

    matched_fragments = []
    matched_fragments_rt = []
    matched_fragments_i = []
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
                        break
    return matched_fragments, matched_fragments_rt, matched_fragments_i

def find_best_peak_group_based_on_reference_sample(d, id, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION):
    for tg in d.keys():
        # build a ngram object for the peak group from reference sample
        ref_pg = build_reference_peak_group(tg, d)

        for sample in id:
            if sample != d[tg].reference_sample.name:
                # for each peak group, create a data structure containing all the information above
                pg = build_other_sample_peak_group(tg, sample, d, ref_pg)

                pg_best = find_best_match_pg(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION)

                # write into display_pg
                d[tg]['display_pg'][sample]['rt_left'] = pg_best['rt_left']
                d[tg]['display_pg'][sample]['rt_right'] = pg_best['rt_right']
                for fragment in pg_best['ms2']['rt_list'].keys():
                    d[tg]['display_pg'][sample]['rt_list'][fragment] = pg_best['ms2']['rt_list'][fragment]
                    d[tg]['display_pg'][sample]['i_list'][fragment] = pg_best['ms2']['i_list'][fragment]
                    d[tg]['display_pg'][sample]['i'][fragment] = pg_best['ms2']['i'][fragment]
                    d[tg]['display_pg'][sample]['if_found_peak'][fragment] = pg_best['ms2']['if_found_peak'][fragment]

                d[tg]['display_pg'][sample]['ms1']['rt_list'] = pg_best['ms1']['rt_list']
                d[tg]['display_pg'][sample]['ms1']['i_list'] = pg_best['ms1']['i_list']
                d[tg]['display_pg'][sample]['ms1']['i'] = pg_best['ms1']['i']
                d[tg]['display_pg'][sample]['ms1']['if_found_peak'] = pg_best['ms1']['if_found_peak']
    return d

def get_fragment_intensity_for_peak_group(peaks_rt, peaks_i, rt, MAX_RT_TOLERANCE):
    rt_dif = 100
    i = 0
    if_found_peak = 0
    for rt0, i0, in zip(peaks_rt, peaks_i):
        rt_dif0 = abs(rt0 - rt)
        if rt_dif0 < rt_dif:
            i = i0
            rt_dif = rt_dif0

    if rt_dif < MAX_RT_TOLERANCE:
        if_found_peak = 1
    return i, if_found_peak

def get_peak_group_values(pg, rt, ref_pg, MAX_RT_TOLERANCE):
    pg_best = data_holder.DataHolder()
    pg_best['peak_rt'] = rt
    # use peak width from ref_pg to compute the boundary of this peak group
    pg_best['rt_left'] = rt - 0.5 * (ref_pg['rt_right'] - ref_pg['rt_left'])
    pg_best['rt_right'] = rt + 0.5 * (ref_pg['rt_right'] - ref_pg['rt_left'])

    for fragment in pg[rt]['ms2']['rt_list'].keys():

        # get the rt and i list only for the peak group
        pg_best['ms2']['rt_list'][fragment], pg_best['ms2']['i_list'][fragment] = \
            chrom.get_chrom_range(pg_best['rt_left'], pg_best['rt_right'],
                                  pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['i_list'][fragment])

        pg_best['ms2']['i'][fragment], pg_best['ms2']['if_found_peak'][fragment] = get_fragment_intensity_for_peak_group(\
            pg[rt]['ms2']['peaks_rt'][fragment], pg[rt]['ms2']['peaks_i'][fragment], rt, MAX_RT_TOLERANCE)

    # get the rt and i list only for the peak group
    pg_best['ms1']['rt_list'], pg_best['ms1']['i_list'] = \
        chrom.get_chrom_range(pg_best['rt_left'], pg_best['rt_right'], pg[rt]['ms1']['rt_list'], pg[rt]['ms1']['i_list'])

    pg_best['ms1']['i'], pg_best['ms1']['if_found_peak'] = get_fragment_intensity_for_peak_group(\
            pg[rt]['ms1']['peaks_rt'], pg[rt]['ms1']['peaks_i'], rt, MAX_RT_TOLERANCE)

    return pg_best

def find_top_n_fragment(option, rt, pg):
    # sort fragment based on i, and then select the top n fragment
    fragment_sorted = []

    for fragment, i in sorted(zip(pg[rt]['ms2']['i'].values(), pg[rt]['ms2']['i'].values()), reverse = True):
        fragment_sorted.append(fragment)

    return fragment_sorted[int(option) - 1]

def filter_peak_group_top_fragment(n, pg, MAX_RT_TOLERANCE):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        if_peak_found = 0
        for rt0 in pg[rt]['ms2']['peaks_rt'][fragment]:
            if abs(rt - rt0) < MAX_RT_TOLERANCE:
                if_peak_found = 1
                break
        if if_peak_found == 0:
            del pg[rt]
    return pg

def filter_peak_group_ms1(pg, MAX_RT_TOLERANCE):
    for rt in pg.keys():
        if_peak_found = 0
        for rt0 in pg[rt]['ms1']['peaks_rt']:
            if abs(rt - rt0) < MAX_RT_TOLERANCE:
                if_peak_found = 1
                break
        if if_peak_found == 0:
            del pg[rt]
    return pg

def compute_transition_intensity (rt0, rt_list, i_list):
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

def filter_peak_group_peak_shape(n, pg, PEAK_SHAPE_FOLD_VARIATION):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        peak_intensity_apex = compute_transition_intensity(rt, pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        peak_intensity_left = compute_transition_intensity(pg[rt]['ms2']['rt_left'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        peak_intensity_right = compute_transition_intensity(pg[rt]['ms2']['rt_right'][fragment], pg[rt]['ms2']['rt_list'][fragment], pg[rt]['ms2']['rt_list'][fragment])
        fold_change_left = float(peak_intensity_apex) / float(peak_intensity_left)
        fold_change_right = float(peak_intensity_apex) / float(peak_intensity_right)

        if  fold_change_left > PEAK_SHAPE_FOLD_VARIATION and fold_change_right > PEAK_SHAPE_FOLD_VARIATION:
            pass  # good peak boundary
        else:
            del pg[rt]
    return pg

def most_correlated_peak_group_based_on_fragment_intensity(pg, ref_pg):
    pg_corr = {}
    for rt in pg.keys():
        x = [ref_pg['ms2']['i'][fragment] for fragment in pg[rt]['ms2']['peaks_i'].keys()]
        y = []
        for fragment in pg[rt]['ms2']['peaks_i'].keys():
            y0 = get_intensity_for_closest_rt(pg[rt]['ms2']['peaks_i'][fragment], pg[rt]['ms2']['peaks_rt'][fragment], rt)
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

def filter_peak_group_top_fragment_peak_boundary(n, pg, ref_pg, PEAK_WIDTH_FOLD_VARIATION):
    for rt in pg.keys():
        fragment = find_top_n_fragment(n, rt, pg)
        ref_sample_peak_width = float(ref_pg['rt_right'] - ref_pg['rt_left'])
        peak_width = float(pg[rt]['ms2']['rt_right'][fragment] - pg[rt]['ms2']['rt_left'][fragment])
        if 1.0 / PEAK_WIDTH_FOLD_VARIATION < peak_width / ref_sample_peak_width < PEAK_WIDTH_FOLD_VARIATION:
            pass  # good peak boundary
        else:
            del pg[rt]
    return pg

def find_top_fragment_with_peak(pg, MAX_RT_TOLERANCE):

    top_fragment = {}

    for rt in pg.keys():
        if_peak_found = 0
        top_n_fragment_found = -1
        n = 1
        for fragment in sorted(pg[rt]['ms2']['i'].keys(), reverse = True):
            for rt0 in pg[rt]['ms2']['peaks_rt'][fragment]:
                if abs(rt0 - rt) < MAX_RT_TOLERANCE:
                    if_peak_found = 1
                    break
            if if_peak_found == 1:
                top_n_fragment_found = n
                break
            else:
                n += 1
        top_fragment[rt] = top_n_fragment_found

    return top_fragment



def find_best_match_pg(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):

    if len(pg) == 0:
        return 0
    elif len(pg) == 1:
        pg_best = only_one_pg( pg, ref_pg )
        return pg_best
    else:
        pg_best = find_best_match_pg_rule_A(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
        return pg_best

def only_one_pg(pg, ref_pg):
    return get_peak_group_values( pg, pg.keys( )[0], ref_pg, MAX_RT_TOLERANCE)


def find_best_match_pg_rule_A(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 1 fragment as a peak
    pg2 = filter_peak_group_top_fragment(1, pg, MAX_RT_TOLERANCE)
    if len(pg2) == 1:
        pg_best = only_one_pg(pg2, ref_pg)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_B(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_B(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best


def find_best_match_pg_rule_B(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 2 fragment as a peak
    pg2 = filter_peak_group_top_fragment(2, pg, MAX_RT_TOLERANCE)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_C(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_C(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best

def find_best_match_pg_rule_C(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without MS1 as a peak
    pg2 = filter_peak_group_ms1(pg, MAX_RT_TOLERANCE)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_D(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_D(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best


def find_best_match_pg_rule_D(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 1 fragment showing good peak boundary
    pg2 = filter_peak_group_top_fragment_peak_boundary(1, pg, ref_pg, PEAK_WIDTH_FOLD_VARIATION)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_E(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_E(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best

def find_best_match_pg_rule_E(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 2 fragment showing good peak boundary
    pg2 = filter_peak_group_top_fragment_peak_boundary(2, pg, ref_pg, PEAK_WIDTH_FOLD_VARIATION)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_F(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_F(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best

def find_best_match_pg_rule_F(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 1 fragment showing good peak shape
    pg2 = filter_peak_group_peak_shape(1, pg, PEAK_SHAPE_FOLD_VARIATION)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_G(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_G(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best

def find_best_match_pg_rule_G(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # filter out peak groups without top 2 fragment showing good peak shape
    pg2 = filter_peak_group_peak_shape(2, pg, PEAK_SHAPE_FOLD_VARIATION)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        pg_best = find_best_match_pg_rule_H(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    elif len(pg2) > 1:
        pg_best = find_best_match_pg_rule_H(pg2, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION)
    return pg_best

def find_best_match_pg_rule_H(pg, ref_pg, MAX_RT_TOLERANCE, PEAK_WIDTH_FOLD_VARIATION, PEAK_SHAPE_FOLD_VARIATION):
    # select the peak group with highest correlation to the reference peak group in terms of fragment intensity
    pg2 = filter_peak_group_peak_shape(2, pg, PEAK_SHAPE_FOLD_VARIATION)
    if len(pg2) == 1:
        pg_best = get_peak_group_values(pg2, pg2.keys()[0], ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) == 0:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(pg, ref_pg)
        pg_best = get_peak_group_values(pg, best_pg_rt, ref_pg, MAX_RT_TOLERANCE)
    elif len(pg2) > 1:
        best_pg_rt = most_correlated_peak_group_based_on_fragment_intensity(pg2, ref_pg)
        pg_best = get_peak_group_values(pg2, best_pg_rt, ref_pg, MAX_RT_TOLERANCE)
    return pg_best


def build_other_sample_peak_group(tg, sample, d, ref_pg):

    pg = data_holder.DataHolder()
    for peak_rt in d[tg]['peak_groups'][sample].keys():
        pg[peak_rt]['ms1']['rt_list'] = d[tg]['precursor']['rt_list'][sample]
        pg[peak_rt]['ms1']['i_list'] = d[tg]['precursor']['i_list'][sample]
        pg[peak_rt]['ms1']['peaks_rt'] = d[tg]['precursor']['peaks_rt'][sample]
        pg[peak_rt]['ms1']['peaks_i'] = d[tg]['precursor']['peaks_i'][sample]

        for fragment in ref_pg['rt_list'].keys():
            pg[peak_rt]['ms2']['rt_list'][fragment] = d[tg]['fragments'][fragment]['rt_list'][sample]
            pg[peak_rt]['ms2']['i_list'][fragment] = d[tg]['fragments'][fragment]['i_list'][sample]
            pg[peak_rt]['ms2']['peaks_rt'][fragment] = d[tg]['fragments'][fragment]['peaks_rt'][sample]
            pg[peak_rt]['ms2']['peaks_i'][fragment] = d[tg]['fragments'][fragment]['peaks_i'][sample]
            pg[peak_rt]['ms2']['rt_left'][fragment] = d[tg]['fragments'][fragment]['peaks_i'][sample]
            pg[peak_rt]['ms2']['i'][fragment] = d[tg]['peak_groups'][sample][peak_rt]['i'][fragment]

    return pg

def build_reference_peak_group(tg, d):
    ref_pg = data_holder.DataHolder()
    ref_sample = d[tg]['reference_sample']['name']
    ref_pg['peak_rt'] = d[tg]['reference_sample']['peak_rt_found']
    ref_pg['rt_left'] = d[tg]['reference_sample']['rt_left']
    ref_pg['rt_right'] = d[tg]['reference_sample']['rt_right']

    ref_pg['ms2']['rt_list'] = d[tg]['display_pg'][ref_sample]['rt_list']
    ref_pg['ms2']['i_list'] = d[tg]['display_pg'][ref_sample]['i_list']
    ref_pg['ms2']['i'] = d[tg]['display_pg'][ref_sample]['i_list']

    ref_pg['ms1']['rt_list'] = d[tg]['display_pg']['precursor']['rt_list']
    ref_pg['ms1']['i_list'] = d[tg]['display_pg']['precursor']['i_list']
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


def find_peak_groups(chrom_data, sample_id, peptide_data):

    peak_groups = {}

    for tg in chrom_data.keys():

        for sample in sample_id:

            all_rt = find_all_rt_values(chrom_data, tg, sample)

            for rt in all_rt:

                peak_groups[rt] = Peak_group(all_rt, chrom_data, peptide_data, tg, sample, rt)

                if len(peak_groups_fragment[rt]) >= MIN_FRAGMENTS: # at least 4 transitions
                    d[tg]['peak_groups'][sample][rt]['fragments'] = peak_groups_fragment[rt]
                    for fragment, i in zip(peak_groups_fragment[rt] , peak_groups_i[rt]):
                        d[tg]['peak_groups'][sample][rt]['i'][fragment] = i
                    for fragment in peak_groups_fragment[rt]:
                        rt_list = map(float, peaks.rt_three_values_to_full_list_string(d[tg]['fragments'][fragment]['rt'][sample]).split(','))
                        i_list = map(float, d[tg]['fragments'][fragment]['i'][sample].split(','))
                        rt_left, rt_right = chrom.get_peak_boundary(rt_list, i_list, rt)
                        d[tg]['peak_groups'][sample][rt]['rt_left'][fragment] = rt_left
                        d[tg]['peak_groups'][sample][rt]['rt_right'][fragment] = rt_right

    return peak_groups


def find_rt_for_reference_sample(d, tg, reference_sample, reference_rt, MAX_RT_TOLERANCE):
    # maybe multiple peak groups are found for the reference sample, here, find the peak rt closest to the rt found by openswath
    good_fragments = []
    rt_dif = MAX_RT_TOLERANCE
    rt_found = -1.0
    for rt in d[tg]['peak_groups'][reference_sample].keys():
        if abs(rt - reference_rt) < rt_dif:
            good_fragments = d[tg]['peak_groups'][reference_sample][rt]['fragments']
            rt_dif = abs(rt - reference_rt)
            rt_found = rt
    return good_fragments, rt_dif, rt_found


def refine_fragments_based_on_reference_sample(d, MAX_RT_TOLERANCE, MIN_FRAGMENTS):

    for tg in d.keys():

        reference_sample = d[tg]['reference_sample']['name']
        reference_rt = d[tg]['reference_sample']['rt']

        good_fragments, rt_dif, rt_found = find_rt_for_reference_sample(d, tg, reference_sample, reference_rt, MAX_RT_TOLERANCE)

        # remove rt for other peak groups, only keep the picked peak group rt
        for rt in d[tg]['peak_groups'][reference_sample].keys():
            if rt != rt_found:
                del d[tg]['peak_groups'][reference_sample][rt]

        # delete data for non-selected fragments
        for fragment in d[tg]['fragments'].keys():
            if fragment in good_fragments:
                pass
            else:
                del d[tg]['fragments'][fragment]

    # remove tg with < 4 good fragments in the reference sample
    for tg in d.keys():
        if len(d[tg]['fragments'].keys()) < MIN_FRAGMENTS:
            del d[tg]
    return d

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