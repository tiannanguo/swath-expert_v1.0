__author__ = 'Tiannan Guo, ETH Zurich 2015'
#this code parses chrom.mzML into chrom.txt file

#compress rt_list to three values: start, end, gap
import sys
import csv
import base64
import struct
import gzip
from collections import defaultdict

class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

# RT_TOLERANCE = 240.0 # 180.0  #3min, -3min ~ 3min. Do not futher select RT range. Use all from OpenswathChromatogramExtraction

# in_file1 = "2R_water.txt"
# in_file2 = "gold1_10.chrom.mzML.gz"
# in_file3 = "AQUA374_step4.csv" # "AQUA374_ms1.csv"  #AQUA374_step4.csv
# out_file = "gold1_10.chrom.txt.gz"
# id_mapping_file = "goldenSets90.txt"

#
in_file1 = sys.argv[1]  #'pcv_reference_peak_groups.txt'
in_file2 = sys.argv[2]  #'PCV1_1.chrom.mzML.gz'
in_file3 = sys.argv[3]  #'140306PC-DDA-92files_step6_pcv.csv'
out_file = sys.argv[4]  #'PCV1_1.chrom.txt.gz'
id_mapping_file = sys.argv[5]  #'goldenSets90.txt'

def decode_rt_list(rt_list):
    decoded = base64.decodestring(rt_list)
    tmp_size = len(decoded)/8
    unpack_format1 = ">%dQ" % tmp_size
    rt_list2 = []

    for tmp in struct.unpack(unpack_format1,decoded):
        tmp_i = struct.pack(">Q", tmp)
        tmp_f = struct.unpack("d", tmp_i)[0]
        rt_list2.append(round(float(tmp_f), 2))

    return rt_list2

def decode_i_list(i_list):

    decoded = base64.decodestring(i_list)
    tmp_size = len(decoded)/4
    unpack_format1 = ">%dL" % tmp_size
    i_list2 = []

    for tmp in struct.unpack(unpack_format1, decoded):
        tmp_i = struct.pack(">I", tmp)
        tmp_f = struct.unpack("f", tmp_i)[0]
        i_list2.append(int(tmp_f))

    return i_list2

def select_rt_range(rt_list, i_list, best_rt, rt_range):

    rt_list2 = [rt for rt in rt_list if best_rt - rt_range < rt < best_rt + rt_range]
    i_list2 = [i for rt, i in zip(rt_list, i_list) if (rt < best_rt + rt_range and rt > best_rt - rt_range)]

    return rt_list2, i_list2

def reduce_rt_list_size(rt_list):

    rt_dif = []

    for rt0, rt1 in zip(rt_list, rt_list[1:]):
        rt_dif.append(float(rt1 - rt0))

    if len(rt_dif) > 0:

        rt_gap = sum(rt_dif) / float(len(rt_dif))
        rt_start = rt_list[0]
        rt_end = rt_list[-1]

        return rt_start, rt_end, rt_gap
    else:
        return 0, 0, 0

def read_reference_peak_group(id_mapping):

    best_rt = {}
    best_sample = {}
    best_score = {}

    with open(in_file1) as CSV_FILE1:

        reader = csv.DictReader(CSV_FILE1, delimiter="\t")

        for row in reader:

            best_rt[row['fragment_name']] = float(row['reference_rt'])
            if bool(id_mapping):
                best_sample[row['fragment_name']] = id_mapping[row['reference_sample'].lower()]
            else:
                best_sample[row['fragment_name']] = row['reference_sample']
            best_score[row['fragment_name']] = float(row['reference_score'])

    return best_rt, best_sample, best_score

def read_library_file():

    library_data = ngram()

    with open(in_file3, 'rb') as CSV_FILE3:
        reader = csv.DictReader(CSV_FILE3, delimiter="\t")

        for row in reader:
            #for MS2/fragments
            library_data[row['transition_name']]['protein'] = row['ProteinName']
            library_data[row['transition_name']]['transition_group_id'] = row['transition_group_id']
            library_data[row['transition_name']]['rt'] = row['Tr_recalibrated']
            library_data[row['transition_name']]['PrecursorMz'] = row['PrecursorMz']
            library_data[row['transition_name']]['ProductMz'] = row['ProductMz']
            library_data[row['transition_name']]['LibraryIntensity'] = row['LibraryIntensity']

            #for MS1/precursors
            library_data[row['transition_group_id']]['protein'] = row['ProteinName']
            library_data[row['transition_group_id']]['transition_group_id'] = row['transition_group_id']
            library_data[row['transition_group_id']]['rt'] = row['Tr_recalibrated']
            library_data[row['transition_group_id']]['PrecursorMz'] = row['PrecursorMz']
            library_data[row['transition_group_id']]['ProductMz'] = row['PrecursorMz']
            library_data[row['transition_group_id']]['LibraryIntensity'] = -1.0   #no library intensity for MS1

    return library_data

def read_id_mapping_file():
    id_mapping = {}

    with open(id_mapping_file, 'rb') as MAPPING_FILE:
        reader = csv.reader(MAPPING_FILE, delimiter="\t")
        for row in reader:
            id_mapping[row[0].lower()] = row[1]

    return id_mapping

def main():

    #check whether there is a id mapping file
    id_mapping = {}
    if bool(id_mapping_file):
        id_mapping = read_id_mapping_file()

    #read best peak group information
    best_rt, best_sample, best_score = read_reference_peak_group(id_mapping)

    #read library file, get peptide ~ protein ~ transition information
    library_data = read_library_file()

    #read chrom.mzML.gz file

    IN_FILE2 = gzip.open(in_file2, 'rb')
    tree = ET.ElementTree(file=IN_FILE2)

    transition_name = ''
    rt_list = ''
    i_list = ''

    #the mzML file contains HUPO_TAG for each tag
    HUPO_TAG = '{http://psi.hupo.org/ms/mzml}'
    chromatogram_tag = ''.join([HUPO_TAG, 'chromatogram'])
    binaryDataArrayList_tag = ''.join([HUPO_TAG, 'binaryDataArrayList'])
    binaryDataArray_tag = ''.join([HUPO_TAG, 'binaryDataArray'])
    binaryParam_tag = ''.join([HUPO_TAG, 'cvParam'])
    binary_tag = ''.join([HUPO_TAG, 'binary'])

    with gzip.open (out_file, 'wb') as OUT_FILE:

        writer = csv.writer(OUT_FILE, delimiter="\t")
        writer.writerow(('transition_name', 'rt_list', 'i_list',
                        'transition_group_id',
                        'best_rt', 'best_sample', 'best_score',
                        'protein', 'irt', 'precursor_mz', 'transition_mz', 'transition_i'))

        for chrom in tree.iter(chromatogram_tag):

            transition_name = chrom.attrib.get('id')

            # some tg is not in best_peak_group therefore remove them
            if library_data[transition_name]['transition_group_id'] in best_rt.keys() :
                # print transition_name
                for chrom in chrom.iter(binaryDataArrayList_tag):
                    for binary_array_data_list in chrom.iter(binaryDataArray_tag):
                        for binary_array_data in binary_array_data_list.iter():
                            if_rt_list = 0
                            if_i_list = 0
                            for cvParam in binary_array_data.iter(binaryParam_tag):
                                if cvParam.attrib.get('name') == 'time array':
                                    if_rt_list = 1
                                elif cvParam.attrib.get('name') == 'intensity array':
                                    if_i_list = 1

                            for binary in binary_array_data.iter(binary_tag):
                                if if_rt_list == 1:
                                    rt_list = binary.text
                                    if_rt_list = 0
                                elif if_i_list == 1:
                                    i_list = binary.text
                                    if_i_list = 0

                                    #decode rt and int
                                    rt_list2 = decode_rt_list(rt_list)
                                    i_list2 = decode_i_list(i_list)

                                    #clear RAM
                                    rt_list = ''
                                    i_list = ''

                                    #select rt range
                                    # rt_list3, i_list3 = select_rt_range(rt_list2, i_list2, best_rt[library_data[transition_name]['transition_group_id']], RT_TOLERANCE)
                                    # do not select rt range here.
                                    rt_list3 = rt_list2
                                    i_list3 = i_list2
                                    # print "rt is %s" %rt_list3
                                    rt_list4 = reduce_rt_list_size(rt_list3)
                                    #clear RAM
                                    rt_list2 = ''
                                    i_list2 = ''
                                    rt_list3 = ''

                                    writer.writerow((transition_name, ','.join(str(p) for p in rt_list4), ','.join(str(p) for p in i_list3),
                                                     library_data[transition_name]['transition_group_id'],
                                                     best_rt[library_data[transition_name]['transition_group_id']], best_sample[library_data[transition_name]['transition_group_id']], best_score[library_data[transition_name]['transition_group_id']],
                                                     library_data[transition_name]['protein'], library_data[transition_name]['rt'], library_data[transition_name]['PrecursorMz'],
                                                     library_data[transition_name]['ProductMz'], library_data[transition_name]['LibraryIntensity']))
                                    #clear RAM
                                    transition_name = ''
                                    rt_list3 = ''
                                    i_list3 = ''

main()