__author__ = 'guot'

import csv
import gzip
import data_holder


def read_sample_replicate_info(sample_replicates_info_file):

    unqiue_sample = {}

    with open(sample_replicates_info_file, 'r') as f:
        line_num = 1
        for line in f:
            unqiue_sample[line_num] = line.split("\t")
            line_num += 1

    return unqiue_sample


def read_id_file(id_mapping_file):
    sample_id = []
    id_mapping = {}
    with open(id_mapping_file) as i:
        reader = csv.reader(i, delimiter="\t")
        for row in reader:
            sample_id.append(row[1])
            id_mapping[row[0].lower()] = row[1]
    return sample_id, id_mapping


def read_com_chrom_file(chrom_file, sample_id, normalization_factors):

    ref_sample_data = {}
    chrom_data = data_holder.NestedDict()
    peptide_data = data_holder.NestedDict()

    # sometimes, in the golden standard data set, multiple "best_sample"s are found. The best in water may not be the best in human
    # in this case, check the input file and write to a new file
    # select the best sample with lowest m_score
    chrom_file2 = use_one_best_sample(chrom_file)

    with gzip.open(chrom_file2, 'rb') as i:
        r = csv.DictReader(i, delimiter='\t')
        for row in r:
            fragment = row['transition_name']
            tg = row['transition_group_id']
            peak_rt = float(row['best_rt'])
            ref_sample_name = row['best_sample']
            ref_sample_score = float(row['best_score'])

            ref_sample_data[tg] = data_holder.ReferenceSample(
                ref_sample_name, ref_sample_score, peak_rt)

            peptide_data[tg]['ms1']['preMz'] = float(row['precursor_mz'])
            peptide_data[tg]['ms1']['protein'] = row['protein']
            peptide_data[tg]['ms1']['irt'] = float(row['irt'])

            peptide_data[tg]['ms2'][fragment] = float(row['transition_mz'])

            for k in sample_id:

                # print k

                if k == "pcf112":
                    pass

                rt_list_three_values_csv = row[k + '_rt']
                i_list_csv = row[k + '_i']

                if rt_list_three_values_csv == "NA" or rt_list_three_values_csv == "0,0,0" or i_list_csv == "NA" or i_list_csv == 0:
                    i_list_csv = "NA"
                    chrom_data[tg][k][fragment] = "NA"

                else:

                    i_list_csv = apply_normalization_based_on_tic(i_list_csv, normalization_factors, k)

                    ###########
                    # if k == 'BR94': # and fragment.startswith('1191_'):
                    #     pass

                    # print k, fragment

                    chrom_data[tg][k][fragment] = data_holder.Chromatogram(
                        rt_list_three_values_csv, i_list_csv)

    # check if reference sample is good


    return ref_sample_data, chrom_data, peptide_data


def apply_normalization_based_on_tic(i_list_csv, normalization_factors, k):

    if i_list_csv != 'NA':
        i_list = map(float, i_list_csv.split(','))
        norm_factor = compute_norm_factor(k, normalization_factors)

        if len(i_list) > 1:
            i_list2 = [round(x * norm_factor, 1) for x in i_list]
            i_list2_str = map(str, i_list2)
            i_list_csv2 = ','.join(i_list2_str)
            return i_list_csv2
        else:
            return i_list_csv
    else:
        return i_list_csv


def compute_norm_factor(k, normalization_factors):
    max_i = max(normalization_factors.values())
    if normalization_factors[k] <= 0:
        print 'error: sample %s has wrong tic value' % k
    norm_factor = float(max_i) / float(normalization_factors[k])
    return norm_factor


def get_tg_list(chrom_file):
    tg_list = {}
    with gzip.open(chrom_file, 'rb') as i:
        r = csv.DictReader(i, delimiter="\t")
        for row in r:
            tg_list[row['transition_group_id']] = 1
    return tg_list.keys()


def get_best_sample_for_each_tg(chrom_file, tg_list):

    best_sample = {}
    best_score = {}
    best_rt = {}

    for tg in tg_list:
        best_sample[tg] = ''
        best_score[tg] = 1.0
        best_rt[tg] = -1.0

        with gzip.open(chrom_file, 'rb') as i:

            r = csv.DictReader(i, delimiter="\t")

            for row in r:
                tg2 = row['transition_group_id']
                if tg2 == tg and float(row['best_score']) < float(best_score[tg]):
                    best_sample[tg] = row['best_sample']
                    best_rt[tg] = row['best_rt']
                    best_score[tg] = float(row['best_score'])

    return best_sample, best_score, best_rt


def use_one_best_sample(chrom_file):

    chrom_file2 = chrom_file.replace('.txt.gz', '_2.txt.gz')

    # get all tg
    tg_list = get_tg_list(chrom_file)

    # get the best sample for each tg
    best_sample, best_score, best_rt = get_best_sample_for_each_tg(chrom_file, tg_list)

    with gzip.open(chrom_file, 'rb') as i, gzip.open(chrom_file2, 'wb') as o:
        r = csv.DictReader(i, delimiter="\t")
        fieldnames1 = r.fieldnames
        w = csv.DictWriter(o, delimiter="\t", fieldnames=fieldnames1)
        w.writeheader()
        # write data
        for row in r:
            tg = row['transition_group_id']
            row['best_sample'] = best_sample[tg]
            row['best_rt'] = best_rt[tg]
            row['best_score'] = best_score[tg]

            # sometimes the rt or i of a sample can be empty, fill with NA
            for col_name in fieldnames1:
                if not row[col_name]:
                    row[col_name] = "NA"
            w.writerow(row)

    return chrom_file2


def read_tic_normalization_file(tic_normalization_file):

    norm_factor = {}

    with open(tic_normalization_file, 'rb') as i:
        r = csv.DictReader(i, delimiter="\t")
        for row in r:
            norm_factor[row['sample']] = row['tic']

    return norm_factor
