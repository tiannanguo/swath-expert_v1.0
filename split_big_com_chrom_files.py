import gzip
import csv
import sys
from collections import defaultdict

class Nested_dict(defaultdict):
    def __init__(self):
        super(Nested_dict, self).__init__(Nested_dict)

    def __deepcopy__(self):
        return self

in_file = sys.argv[1]
# in_file = 'com_chrom_31.txt.gz'

num_tg_per_file = 30

def get_num_fractions(n):

    print 'totally ' + str(n) + " tg IDs"
    m = n / num_tg_per_file
    res = n - m * num_tg_per_file
    k = 0
    if res == 0:
        k = m
    else:
        k = m + 1
    print 'divide into ' + str(k) + 'fractions'
    return k

def read_data(in_file):

    data = defaultdict(list)

    tg_list = []

    title = ''

    with gzip.open(in_file, 'rb') as IN:
        title = IN.readline()
        for row in IN.readlines():
            tg = row.split("\t")[1]
            data[tg].append(row)
            tg_list.append(tg)

    tg_list = list(set(tg_list))

    num_fractions = get_num_fractions(len(tg_list))

    tg_dict = get_tg_dict(tg_list, num_fractions)

    return title, data, num_fractions, tg_dict

def get_tg_dict(tg_list, num_fractions):

    tg_dict = defaultdict(list)
    last_index = -1

    for i in range(num_fractions - 1):
        start_idx = num_tg_per_file * i
        end_idx = start_idx + num_tg_per_file - 1
        last_index = end_idx
        for j in range(start_idx, end_idx + 1):
            tg_dict[i].append(tg_list[j])

    for j in range(last_index + 1, len(tg_list)):
        tg_dict[num_fractions - 1].append(tg_list[j])

    return tg_dict

def write_data(in_file, data, num_fractions, tg_dict, title):
    for i in range(num_fractions):
        out_file = in_file.replace('.txt.gz', '_' + str(i) + '.txt.gz')
        with gzip.open(out_file, 'wb') as OUT:
            OUT.write(title)
            for tg in tg_dict[i]:
                for d in data[tg]:
                    OUT.write(d)


def main():
    title, data, num_fractions, tg_dict = read_data(in_file)

    write_data(in_file, data, num_fractions, tg_dict, title)

main()