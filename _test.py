__author__ = 'Tiannan Guo, ETH Zurich 2015'
import numpy as np
from math import pi, log
import pylab
from scipy import fft, ifft
from scipy.optimize import curve_fit
from collections import defaultdict



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
    max_peaks_x = []
    max_peaks_y = []

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
    mx = -np.Inf

    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead],
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x

        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks_x.append(mxpos)
                max_peaks_y.append(mx)
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue

    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks_x.pop(0)
            max_peaks_y.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass

    return max_peaks_x, max_peaks_y
    #example of max_peaks: [[1229.1099999999999, 329439.0], [1338.3599999999999, 15618.0], [1389.5699999999999, 42400.0], [1468.0899999999999, 12652.0], [1550.02, 5779.0], [1635.3699999999999, 66082.0], [1833.3800000000001, 111606.0], [1928.97, 3993.0]]



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



def rt_three_values_to_full_list_string(rt):
    rt2 = map(float, rt.split(","))
    rt3 = list (np.arange(rt2[0], rt2[1], rt2[2]))
    #soemtimes the last data pint is lost due to small difference <1 second
    if (rt2[1] - rt3[-1] > 1):
        rt3.append(rt2[1])
    rt4 = [round(s,2) for s in rt3]
    rt5 = map(str, rt4)
    return ",".join(rt5)


def main():
    rt = '1453.01,1932.7,3.42635714286'

    i = '0,0,10,0,30,0,0,0,0,0,0,10,0,20,0,10,0,10,10,10,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,10,0,0,0,0,20,0,10,0,0,20,10,0,0,20,90,61,90,40,10,10,0,0,0,10,0,0,0,0,0,70,181,449,1182,3727,3963,2757,991,182,20,0,20,10,0,0,10,0,30,0,0,0,0,20,20,10,0,0,0,0,0,0,0,10,10,0,0,30,20,10,0,0,0,0,0,0,0,0,0,0,10,0,10,0,0,10,0,10,10,0,0,0,0,0,0,0,0,0,0,0,10,0,20,0,30,0,10'

    rt = map(float, rt_three_values_to_full_list_string(rt).split(','))
    i = map(float, i.split(','))

    print rt

    rt, i = extend_rt_border(rt, i)
    peaks = peakdetect(i, rt, 9.0, 0.3)

    print peaks


main()
