__author__ = 'Tiannan Guo, ETH Zurich 2015'
import gzip
import csv


png_id = 205

def get_file_num(png_id):
    target_file = ''
    if_found = 0
    for n in range(1,33):
        chrom_file = 'com_chrom_%s.txt.gz' % str(n)
        with gzip.open(chrom_file, 'rb') as i:
            r = csv.DictReader(i, delimiter='\t')
            for row in r:
                if row['transition_group_id'].split('_')[0] == str(png_id):
                    target_file = chrom_file
                    print 'found file %s' % target_file
                    if_found = 1
                    break
        if if_found == 1:
            break
    return target_file




in_file = get_file_num(png_id)

# in_file = 'com_chrom_%s.txt.gz' % file_num

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




