__author__ = 'Tiannan Guo, ETH Zurich 2015'

import gzip
import csv
import chrom
import peaks
from collections import defaultdict

PLOT_LINE_WIDTH = 2.0  # the thickness of chrom curve

class ngram(defaultdict):
    def __init__(self):
        super(ngram, self).__init__(ngram)

def write_r_code_for_all_samples(d, sample_id):

    all_r_code_samples = []

    for tg in d.keys():

        num_transition = len(d[tg]['fragments'].keys())
        max_intensity_ms1 = get_max_ms1_intensity_in_all_samples(d, tg)
        max_intensity_ms2 = get_max_ms2_intensity_in_all_samples(d, tg)

        r_code_create_png_file = create_png_file(tg, 600, 1000, num_transition)
        all_r_code_samples.append(r_code_create_png_file)

        for sample in sample_id:

            r_code_samples_par = write_sample_par(sample)

            r_code_sample_ms1 = write_sample_ms1(d, sample, tg)

            r_code_sample_ms2, transition_color_code_mapping = write_sample_ms2(d, sample, tg)

            r_code_this_sample = r_code_samples_par + r_code_sample_ms1 + r_code_sample_ms2

            # next step: write out quant table.
            # another expert rules to check the quality of the displayed peak group for quantifiation
            # debug the codes in peak_groups.py, chrom.py, peaks.py, io_swath.py, etc.
            all_r_code_samples.append(r_code_this_sample)

        r_code_close_png_file = write_r_code_close_png_file(tg, num_transition, max_intensity_ms1, max_intensity_ms2, d[tg]['display_pg'][sample_id[0]]['rt_list'].keys()) #contains legend

        all_r_code_samples.append(r_code_close_png_file)

    return all_r_code_samples

def get_max_ms2_intensity_in_all_samples(d, tg):
    i = []
    for sample in d[tg]['display_pg'].keys():
        i0 = max (d[tg]['display_pg'][sample]['i'].values())
        i.append(i0)
    return max(i)

def get_max_ms1_intensity_in_all_samples(d, tg):
    i = []
    for sample in d[tg]['display_pg'].keys():
        i.append(d[tg]['display_pg'][sample]['ms1']['i'])
    return max(i)

def write_r_code_close_png_file(this_tg, num_transitions, max_intensity_ms1, max_intensity_ms2, transition_list):

    r_code_close_file = []

    r_code_close_file.append(
        '''mtext("Peptide precursor id = %s, number of transitions = %i, max MS1 intensity = %.1f, max MS2 intensity = %.1f\n", NORTH<-3, outer= TRUE, cex = 4)\n''' % (
            this_tg, num_transitions, max_intensity_ms1, max_intensity_ms2))


    #add legend
    transition_list_quoted_csv_string = '"'
    for this_transition in transition_list[:-1]:
        transition_list_quoted_csv_string += (this_transition + '", "')
    transition_list_quoted_csv_string += (transition_list[-1] + '"')

    color_code_csv_string = ''
    for n in range(len(transition_list)):
        color_code_csv_string += (str(n + 1) + ', ')
    color_code_csv_string += str(len(transition_list) + 1)

    r_code_close_file.append('''add_legend("bottom", legend = c(%s), pch = 20, pt.cex = 5, cex = 2, horiz=TRUE, col = c(%s))\n'''
                             % (transition_list_quoted_csv_string, color_code_csv_string))
    r_code_close_file.append("dev.off()\n")

    return r_code_close_file

def write_add_legend_function():
    #a function to write legend
    r_code = """add_legend <- function(...) {\n"""
    r_code += """    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)\n"""
    r_code += """    on.exit(par(opar))\n"""
    r_code += """    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')\n"""
    r_code += """    legend(...)\n"""
    r_code += """}\n\n"""

    return r_code


def create_png_file (tg_id, width, height, num_colors):
    code = '''# %s \n''' % tg_id
    code += '''png("%s.png", width = %i, height = %i)\n''' % \
            (tg_id.replace('(', '_').replace(')', '_').replace(':', "_"), width, height)  # change 55_AAAGEFADDPC(UniMod:4)SSVK_2  to 55_AAAGEFADDPC_UniMod_4_SSVK_2

    code += '''par(mfrow=c(2,60), oma = c(5,5,15,5))\n'''
    code += '''palette(rainbow(%d))\n''' % num_colors

    return code


def write_sample_par (sample):
    code = '''#sample %s\n''' % sample
    code += '''par(mar=c(5,1,5,1))\n'''
    return code

def write_sample_ms1(d, sample, tg):

    r_code = '''#MS1 chrom\n'''
    r_code += '''#MS1 id is %s\n''' % tg
    r_code += '''rt = c(%s)\n''' % ','.join(map(str, d[tg]['display_pg'][sample]['ms1']['rt_list']))
    r_code += '''int = c(%s)\n''' % ','.join(map(str, [x * (-1) for x in d[tg]['display_pg'][sample]['ms1']['i_list']]))

    r_code += '''plot(rt, int, type = "l", xlim = c (%.1f, %.1f), ylim = c(-100, 100), ''' \
                         '''lwd = %.1f, xlab = "rt (s)", ylab = "intensity (''' % (
                             d[tg]['display_pg'][sample]['rt_left'],
                             d[tg]['display_pg'][sample]['rt_right'],
                             PLOT_LINE_WIDTH) + '%' + ''')", yaxt = "n", cex.axis = 1.5, frame.plot=FALSE)\n'''
    r_code += '''rm (rt, int)\n'''

    # write text
    # the name of cell line is most case too long for the plot, therefore, wrap to new line for every 5 letters
    this_cell_label = sample
    r_code += '''title("%s", cex.main = 2)\n''' % this_cell_label

    return r_code


def write_sample_ms2(d, sample, tg):

    r_code = ''

    transition_color_code = 1.0
    transition_color_code_mapping = {}

    for fragment in d[tg]['display_pg'][sample]['rt_list'].keys():
        r_code += '''#MS2 chrom\n'''
        r_code += '''#MS2 id is %s\n''' % fragment
        r_code += '''rt = c(%s)\n''' % ','.join(map(str, d[tg]['display_pg'][sample]['rt_list'][fragment]))
        r_code += '''int = c(%s)\n''' % ','.join(map(str, d[tg]['display_pg'][sample]['i_list'][fragment]))
        r_code += '''lines(rt, int, type = "l", col = %d, lwd = %.1f)\n''' % (
            transition_color_code, PLOT_LINE_WIDTH)
        r_code += '''rm (rt, int)\n'''

        transition_color_code_mapping[transition_color_code] = fragment
        transition_color_code += 1

    return r_code, transition_color_code_mapping

def multiple_line_text(label):
    label2 = []
    for i in range(len(label)):
        if i%5 == 0 and i >0 :
            label2.append("\n" + label[i])
        else:
            label2.append(label[i])
    return ''.join(label2)

