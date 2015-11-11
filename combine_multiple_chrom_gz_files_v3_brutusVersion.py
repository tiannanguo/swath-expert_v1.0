__author__ = 'Tiannan Guo, ETH Zurich 2015'

import glob
import re
import gzip
import csv
import time
import sys
import tiannan_handy_defs as tiannan



lib_part_number_dir = sys.argv[1]   #libPart_1, input only 1
size_of_chunk = sys.argv[2]         #size of each chunk. The last chunk might be smaller
chunk_number_to_process = sys.argv[3]    # divide the libPart_2 into equally sized chunks, process one of them here

start_time = time.time()


#select peptides
selected_transition = {}
lib_part_csv_file_name = "LibPart_" + lib_part_number_dir + ".csv"
num_select = 0
with open(lib_part_csv_file_name, 'rb') as LIB_FILE:
    reader = csv.DictReader (LIB_FILE, delimiter = "\t")
    rows = list(reader)
    chunks = tiannan.chunks (rows, int(size_of_chunk))
    for row in chunks[int(chunk_number_to_process)]:
        selected_transition[row['transition_name']] = 1
        selected_transition[row['transition_group_id']] = 1
        num_select += 1

print 'we selected %s transitions to extract in this run' % num_select

#get chrom.gz files
glob_string = "./libPart_" + lib_part_number_dir +"/*.gz"
print "glob ", glob_string
chrom_gz_file_list = glob.glob(glob_string)
#['./libPart_1\\guot_L130610_003_SW_libPart_1.mzML.chrom.gz',
#['./libPart_1/guot_L130610_003_SW_libPart_1.mzML.chrom.gz',
in_file_pattern = re.compile (r'./(libPart_\d+)/(.*)_libPart_\d+.mzML.chrom.gz')
m = in_file_pattern.search(chrom_gz_file_list[0])
lib_part_number = m.group(1)
out_file_name = lib_part_number + "_fraction" + chunk_number_to_process  + ".allChrom.gz"

with gzip.open(out_file_name, 'wb') as OUT_FILE:

    writer = csv.writer(OUT_FILE, delimiter = "\t")

    transition_info_list ={}
    rt_list = {}
    int_list = {}

    in_file_name_short_list =[]

    for in_file_name in chrom_gz_file_list :
        m = in_file_pattern.search(in_file_name)
        in_file_name_short = m.group(2)
        # print "short file name is ", in_file_name_short  #########
        in_file_name_short_list.append(in_file_name_short)
        with gzip.open(in_file_name, 'rb') as IN_FILE:
            reader = csv.DictReader(IN_FILE, delimiter = "\t")
            for row in reader :
                if row['transition_name'] in selected_transition.keys():
                    transition_info_list[row['transition_name']] = (row['transition_group_id'], row['best_rt'], row['best_sample'], row['best_score'], row['protein'], row['irt'], row['precursor_mz'], row['transition_mz'], row['transition_int'])
                    data_id = row['transition_name'] + '_' + in_file_name_short
                    # print "read gzip file, data_id is ", data_id  #############
                    rt_list[data_id] = row['rt_list']
                    int_list[data_id] = row['int_list']

    #write out_file header
    title_to_print =[]
    title_to_print.append('transition_name')
    title_to_print.append('transition_group_id')
    title_to_print.append('best_rt')
    title_to_print.append('best_sample')
    title_to_print.append('best_score')
    title_to_print.append('protein')
    title_to_print.append('irt')
    title_to_print.append('precursor_mz')
    title_to_print.append('transition_mz')
    title_to_print.append('transition_int')
    for in_file_name in in_file_name_short_list:
        title_to_print.append(in_file_name + "_rt")
        title_to_print.append(in_file_name + "_int")


    writer.writerow(title_to_print)


    #print content
    for transition_name in transition_info_list.keys() :
        # print "transition name is ", transition_name   #############
        #rt_list and int_list for all files
        print_string =[]
        print_string.append(transition_name)
        for info in transition_info_list[transition_name]:
            print_string.append(info)

        for in_file_name in in_file_name_short_list :

            data_id = transition_name + '_' + in_file_name

            # print "data id  2 is ", data_id
            if data_id in rt_list.keys():
                print_string.append(rt_list[data_id])
                # print "found data id  2 is ", data_id
            else:
                print_string.append("NA")
                # print "CANNOT find data id  2 is ", data_id

            if data_id in int_list.keys():
                print_string.append(int_list[data_id])
            else:
                print_string.append("NA")

        writer.writerow(print_string)

print("--- %s seconds ---" % (time.time() - start_time))