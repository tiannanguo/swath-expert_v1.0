import numpy as np
from math import pi, log
import pylab
from scipy import fft, ifft
from scipy.optimize import curve_fit
from collections import defaultdict
import parameters
import chrom

def rt_three_values_to_full_list_string(rt):
    rt2 = map(float, rt.split(","))
    rt3 = list (np.arange(rt2[0], rt2[1], rt2[2]))
    #soemtimes the last data pint is lost due to small difference <1 second
    if (rt2[1] - rt3[-1] > 1):
        rt3.append(rt2[1])
    rt4 = [round(s,2) for s in rt3]
    rt5 = map(str, rt4)
    return ",".join(rt5)

def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    
    if len(y_axis) != len(x_axis):
        raise (ValueError, 
                'Input vectors y_axis and x_axis must have same length')
    
    #needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis

def peakdetect(y_axis, x_axis = None, lookahead = 300, delta=0):
    """
    Converted from/based on a MATLAB script at:
    http://billauer.co.il/peakdet.html

    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively

    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200)
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function

    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value)
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do:
        x, y = zip(*tab)
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false

    # check input data
    x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)


    #perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"

    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf

    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead],
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x

        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]

        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]


    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass

    return [max_peaks, min_peaks]
    #example of max_peaks: [[1229.1099999999999, 329439.0], [1338.3599999999999, 15618.0], [1389.5699999999999, 42400.0], [1468.0899999999999, 12652.0], [1550.02, 5779.0], [1635.3699999999999, 66082.0], [1833.3800000000001, 111606.0], [1928.97, 3993.0]]
    # other example:
    # max  [[1692.8599999999999, 3963.0]]
    # min  [[1473.5699999999999, 0.0]]

def extend_rt_border(rt, i):
    mean_gap = 0
    for j in range(1, len(rt)):
        mean_gap += rt[j] - rt[j-1]
    mean_gap = mean_gap / (len(rt) - 1)

    rt_extended_left = np.arange(rt[0] - 5* mean_gap, rt[0] - mean_gap+1, mean_gap)  #1 means extend it a bit because sometimes one data point is missing
    rt_extended_right = np.arange(rt[-1] + mean_gap, rt[-1] + 5* mean_gap, mean_gap)

    int_extended_left = len(rt_extended_left) * [0.0]
    int_extended_right = len(rt_extended_right) * [0.0]

    rt2 = list(rt_extended_left) + rt + list(rt_extended_right)
    i2 = int_extended_left + i + int_extended_right

    return rt2, i2


def select_peak_rt_from_detected_peaks(peak_max):
    #example of max_peaks: [[1229.1099999999999, 329439.0], [1338.3599999999999, 15618.0], [1389.5699999999999, 42400.0], [1468.0899999999999, 12652.0], [1550.02, 5779.0], [1635.3699999999999, 66082.0], [1833.3800000000001, 111606.0], [1928.97, 3993.0]]
    peaks_rt = []
    peaks_i = []
    for pair in peak_max:
        peaks_rt.append(pair[0])
        peaks_i.append(pair[1])
    return peaks_rt, peaks_i


def find_peaks_rt_i(rt, i):

    rt, i = extend_rt_border(rt, i)

    peak_max, peak_min = peakdetect(i, rt, 9.0, 0.3)  #NOTE: input intensity first, then rt

    peaks_rt, peaks_i = select_peak_rt_from_detected_peaks(peak_max)  #NOTE: this includes this_tg, MS1

    return peaks_rt, peaks_i


def find_peaks(d, sample_id):
    for tg in d.keys():
        # ms1
        for sample in sample_id:
            rt = map(float, rt_three_values_to_full_list_string(d[tg]['precursor']['rt_list'][sample]).split(','))
            i = map(float, d[tg]['precursor']['i_list'][sample].split(','))
            if rt[0] <= 0: #sometimes there is no signal
                d[tg]['precursor']['peaks_rt'][sample] = 'NA'
                d[tg]['precursor']['peaks_i'][sample] = 'NA'
            else:
                peaks_rt, peaks_i = find_peaks_rt_i(rt, i)
                if len(peaks_rt) > 0:
                    d[tg]['precursor']['peaks_rt'][sample] = peaks_rt
                    d[tg]['precursor']['peaks_i'][sample] = peaks_i
                else: #no peak found
                    d[tg]['precursor']['peaks_rt'][sample] = 'NA'
                    d[tg]['precursor']['peaks_i'][sample] = 'NA'
        # ms2
        for fragment in tg.keys():
            for sample in sample_id:
                rt = map(float, rt_three_values_to_full_list_string(d[tg]['fragments'][fragment]['rt_list'][sample]).split(','))
                i = map(float, d[tg]['fragments'][fragment]['i_list'][sample].split(','))
                if rt[0] <= 0: #sometimes there is no signal
                    d[tg]['fragments'][fragment]['peaks_rt'][sample] = 'NA'
                    d[tg]['fragments'][fragment]['peaks_i'][sample] = 'NA'
                else:
                    peaks_rt, peaks_i = find_peaks_rt_i(rt, i)
                    if len(peaks_rt) > 0:
                        d[tg]['fragments'][fragment]['peaks_rt'][sample] = peaks_rt
                        d[tg]['fragments'][fragment]['peaks_i'][sample] = peaks_i
                    else:
                        d[tg]['fragments'][fragment]['peaks_rt'][sample] = 'NA'
                        d[tg]['fragments'][fragment]['peaks_i'][sample] = 'NA'
    return d
