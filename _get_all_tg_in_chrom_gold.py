import csv
import gzip

tg = set()
for i in range(1, 33):
    with gzip.open('com_chrom_%d.txt.gz' % i, 'rb') as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            tg.add(row['transition_group_id'])


with open('_all_tg_gold90.txt', 'wb') as fh:
    for tg_id in tg:
        print tg_id
        fh.write(tg_id + "\n")
