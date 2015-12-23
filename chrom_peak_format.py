__author__ = 'guot'

import numpy as np


def rt_three_values_to_full_list_float(rt_list):
    """
    rt time is stored as three values, start, end, interval
    This function converts it into a list of float rounded two to digits after the decimal
    point, csv string for r code
    """
    rt_list2 = map(float, rt_list.split(","))
    rt_list3 = list(np.arange(rt_list2[0], rt_list2[1], rt_list2[2]))
    # sometimes the last data point is lost due to small difference <1 second
    if (rt_list2[1] - rt_list3[-1] > 1):
        rt_list3.append(rt_list2[1])
    rt_list4 = [round(s, 2) for s in rt_list3]
    return rt_list4


def rt_three_values_to_full_list_string(rt_list):
    rt_list2 = map(float, rt_list.split(","))
    rt_list3 = list(np.arange(rt_list2[0], rt_list2[1], rt_list2[2]))
    # soemtimes the last data pint is lost due to small difference <1 second
    if (rt_list2[1] - rt_list3[-1] > 1):
        rt_list3.append(rt_list2[1])
    rt_list4 = [round(s, 2) for s in rt_list3]
    rt_list5 = map(str, rt_list4)
    return ",".join(rt_list5)


def extend_rt_border(rt_list, int_list):
    """when detect peaks, for most normal cases, the peakdetect.py works well.
    but for some peaks at the border of which one boundary is cut out, this function can not detect
    the peak.  in this case, this function try to repair the peak.
    """
    mean_gap = 0
    num_data = 0
    for i in range(1, len(rt_list)):
        mean_gap += rt_list[i] - rt_list[i - 1]
        num_data += 1
    mean_gap = mean_gap / num_data

    # 1 means extend it a bit because sometimes one data point is missing
    rt_extended_left = np.arange(rt_list[0] - 5 * mean_gap, rt_list[0] - mean_gap + 1, mean_gap)
    rt_extended_right = np.arange(rt_list[-1] + mean_gap, rt_list[-1] + 5 * mean_gap, mean_gap)

    int_extended_left = len(rt_extended_left) * [0.0]
    int_extended_right = len(rt_extended_right) * [0.0]
    #
    # print len(rt_extended_left), len(rt_extended_right), len(int_extended_left), len(int_extended_right)
    #
    # print len(rt_list), len(int_list)

    rt2 = list(rt_extended_left) + rt_list + list(rt_extended_right)
    int2 = int_extended_left + int_list + int_extended_right

    # print len(rt2), len(int2)

    return rt2, int2


def select_peaks_from_detected_peaks(peak_max):
    """
    the results from peakdetect is a list of tuple, this function is to take
    out the rt information of peaks I need
    # example of max_peaks: [[1229.1099999999999, 329439.0],
    # [1338.3599999999999, 15618.0], [1389.5699999999999, 42400.0],
    # [1468.0899999999999, 12652.0], [1550.02, 5779.0], [1635.3699999999999,
    # 66082.0], [1833.3800000000001, 111606.0], [1928.97, 3993.0]]
    """
    this_all_peaks_list = []
    for peak_max_pair in peak_max:
        # print "type of peak_max_rt is ", type(peak_max_pair[0])
        this_all_peaks_list.append(peak_max_pair[0])
    return this_all_peaks_list


def cut_rt_range_for_display_quant(rt_list, int_list, int_list2, rt_left, rt_right):
    rt_list2 = []
    int_list3 = []
    int_list4 = []
    for rt, i, i2 in zip(rt_list, int_list, int_list2):
        if rt_left < rt < rt_right:
            rt_list2.append(rt)
            int_list3.append(i)
            int_list4.append(i2)
    return rt_list2, int_list3, int_list4


def format_and_normalize_intensity(int_list, max_int, MS1_MS2):
    """
    convert the intensity values into normalized in -100 to 100, and then
    save in a csv string for r code
    """
    int_list2 = map(float, int_list.split(","))
    MS1_or_MS2 = 1.0
    if (MS1_MS2 == 'MS1'):
        MS1_or_MS2 = -1.0
    int_list3 = [i / max_int * 100 * MS1_or_MS2 for i in int_list2]
    int_list4 = [round(s, 2) for s in int_list3]
    return int_list2, int_list4
