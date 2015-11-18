__author__ = 'Tiannan Guo, ETH Zurich 2015'
import gzip
import csv

file_num = 17
png_id = 291

in_file = 'com_chrom_%s.txt.gz' % file_num
out_file = 'debug_png_id_%s.txt.gz' % png_id

with gzip.open(in_file, 'rb') as i, gzip.open(out_file, 'wb') as o:
    r = csv.DictReader(i, delimiter="\t")
    fieldnames1 = r.fieldnames
    w = csv.DictWriter(o, delimiter="\t", fieldnames=fieldnames1)
    w.writeheader()
    for row in r:
        tg = row['transition_group_id']
        tg2 = tg.split('_')
        tg_id = tg2[0]
        if str(tg_id) == str(png_id):
            w.writerow(row)




