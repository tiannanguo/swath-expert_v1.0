__author__ = 'Tiannan Guo'
"""
Tiannan Guo, ETH Zurich, 2015
"""
import re
import csv
import sys

input_file_name = sys.argv[1]   #"feature_alignment_requant_matrix_test.tsv"

output_file_name = sys.argv[2]   # "best_peak_group.txt"

with open(input_file_name,'rb') as csvfile, open(output_file_name, 'wb') as out_file:

    reader = csv.DictReader(csvfile, delimiter="\t")
    writer = csv.writer(out_file, delimiter="\t")

    writer.writerow(('fragment_name', 'reference_sample', 'reference_score', 'reference_rt'))

    for row in reader:
        best_score = 1.0
        best_sample = ''
        best_rt = 0.0
        m1 = re.search(r'^(\d+\_.*)\_run0', row['Peptide'])
        if m1:
            for column_name in row.keys():

                m2 = re.search(r'score_(.*?)(_with_dscore.*_\d+\_\d+)', column_name)

                if m2:

                    if row[column_name] != 'NA':
                        sample_id = m2.group(1)
                        rt_column_name = "RT_" + sample_id + m2.group(2)
                        this_score = float(row[column_name])

                        if this_score < best_score :

                            best_score = float(row[column_name])
                            best_sample = sample_id
                            best_rt = float(row[rt_column_name])

            writer.writerow((row["Peptide"][:-5], best_sample, best_score, best_rt))
