import csv
import gzip

tg = {}
for i in range(1,33):
    with gzip.open('com_chrom_' + str(i) + '.txt.gz', 'rb') as i:
        r = csv.DictReader(i, delimiter="\t")

        for row in r:

            tg[row['transition_group_id']] = 1


with open('_all_tg_gold90.txt', 'wb') as OUT:

    for tg_id in tg.keys():
        print tg_id
        OUT.write(tg_id + "\n")