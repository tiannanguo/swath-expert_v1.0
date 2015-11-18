import csv
import gzip

chrom_file = 'debug_png_id_526.txt.gz'

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
            tg = row['transition_group_id']
            if float(row['best_score']) < float(best_score[tg]):
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
            w.writerow(row)

    return chrom_file2

def main():
    use_one_best_sample(chrom_file)

main()

