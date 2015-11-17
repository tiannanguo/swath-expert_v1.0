__author__ = 'Tiannan Guo, ETH Zurich 2015'

import gzip
import csv
import chrom
import peaks
from collections import defaultdict
import parameters

def write_r_code_for_all_samples(display_data, sample_id, out_R_file):

    all_r_code_samples = []

    for tg in display_data.keys():

        num_transition = len(display_data[tg][sample_id[0]]['ms2']['rt_list'].keys())
        max_intensity_ms1 = get_max_ms1_intensity_in_all_samples(display_data, tg)
        max_intensity_ms2 = get_max_ms2_intensity_in_all_samples(display_data, tg)

        r_code_create_png_file = create_png_file(tg, parameters.PNG_FILE_WIDTH, parameters.PNG_FILE_HEIGHT, num_transition)
        all_r_code_samples.append(r_code_create_png_file)

        for sample in sample_id:

            r_code_samples_par = write_sample_par(sample)

            r_code_sample_ms1 = write_sample_ms1(display_data, sample, tg, max_intensity_ms1)

            r_code_sample_ms2, transition_color_code_mapping = write_sample_ms2(display_data, sample, tg, max_intensity_ms2)

            r_code_this_sample = r_code_samples_par + r_code_sample_ms1 + r_code_sample_ms2

            # next step: write out quant table.
            # another expert rules to check the quality of the displayed peak group for quantifiation
            # debug the codes in peak_groups.py, chrom.py, peaks.py, io_swath.py, etc.
            all_r_code_samples.append(r_code_this_sample)

        r_code_close_png_file = write_r_code_close_png_file(
            tg, num_transition, max_intensity_ms1, max_intensity_ms2,
            display_data[tg][sample_id[0]]['ms2']['rt_list'].keys()) #contains legend

        all_r_code_samples.append(r_code_close_png_file)


    #TODO write R codes into the file out_R_file
    with open(out_R_file, 'wb') as out_file:
        for line in all_r_code_samples:
            out_file.write(line)

    return 1

def get_max_ms2_intensity_in_all_samples(display_data, tg):

    i = []

    for sample in display_data[tg].keys():

        all_i = display_data[tg][sample]['ms2']['peak_apex_i'].values()

        if len(all_i) > 0:
            i0 = max(all_i)
            i.append(i0)

    max_i = float(max(i))

    if max_i == 0:
        max_i = 0.1

    return max_i

def get_max_ms1_intensity_in_all_samples(display_data, tg):

    i = []
    for sample in display_data[tg].keys():
        i.append(display_data[tg][sample]['ms1']['peak_apex_i'])

    max_i = float(max(i))

    if max_i == 0:
        max_i = 0.1

    return max_i

def write_r_code_close_png_file(this_tg, num_transitions, max_intensity_ms1, max_intensity_ms2, transition_list):

    r_code_close_file = []

    r_code_close_file.append(
        '''mtext("Peptide precursor id = %s, number of transitions = %i, max MS1 intensity = %.1f, max MS2 intensity = %.1f\n", NORTH<-3, outer= TRUE, cex = %s)\n''' % (
            this_tg, num_transitions, max_intensity_ms1, max_intensity_ms2, parameters.title_font_size))

    # add legend
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

    return ''.join(r_code_close_file)

def write_add_legend_function():
    #a function to write legend
    r_code = """add_legend <- function(...) {\n"""
    r_code += """    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)\n"""
    r_code += """    on.exit(par(opar))\n"""
    r_code += """    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')\n"""
    r_code += """    legend(...)\n"""
    r_code += """}\n\n"""

    return r_code


def create_png_file(tg_id, width, height, num_colors):

    code = '''# %s \n''' % tg_id
    code += '''png("%s.png", width = %i, height = %i)\n''' % \
            (tg_id.replace('(', '_').replace(')', '_').replace(':', "_"), width, height)  # change 55_AAAGEFADDPC(UniMod:4)SSVK_2  to 55_AAAGEFADDPC_UniMod_4_SSVK_2

    code += '''par(mfrow=c(%s,%s), oma = c(5,5,15,5))\n''' % (parameters.figures_num_rows, parameters.figure_num_per_row)
    code += '''palette(rainbow(%d))\n''' % num_colors

    return code


def write_sample_par(sample):
    code = '''#sample %s\n''' % sample
    code += '''par(mar=c(5,1,5,1))\n'''
    return code

def write_sample_ms1(display_data, sample, tg, max_intensity_ms1):

    r_code = '''#MS1 chrom\n'''
    r_code += '''#MS1 id is %s\n''' % tg

    rt_list = ','.join(map(str, display_data[tg][sample]['ms1']['rt_list']))

    r_code += '''rt = c(%s)\n''' % rt_list

    i_list = ','.join(map(str, [round(x * (-100) / max_intensity_ms1, 1) for x in display_data[tg][sample]['ms1']['i_list']]))

    r_code += '''int = c(%s)\n''' % i_list

    r_code += '''plot(rt, int, type = "l", xlim = c (%.1f, %.1f), ylim = c(-100, 100), ''' \
                         '''lwd = %.1f, xlab = "rt (s)", ylab = "intensity (''' % (
                             display_data[tg][sample]['rt_left'],
                             display_data[tg][sample]['rt_right'],
                             parameters.PLOT_LINE_WIDTH) + '%' + ''')", yaxt = "n", cex.axis = 1.5, frame.plot=FALSE)\n'''
    r_code += '''rm (rt, int)\n'''

    # write text
    # the name of cell line is most case too long for the plot, therefore, wrap to new line for every 5 letters
    this_cell_label = sample
    r_code += '''title("%s", cex.main = 2)\n''' % this_cell_label

    return r_code


def write_sample_ms2(display_data, sample, tg, max_intensity_ms2):

    r_code = ''

    transition_color_code = 1.0
    transition_color_code_mapping = {}

    for fragment in display_data[tg][sample]['ms2']['rt_list'].keys():
        r_code += '''#MS2 chrom\n'''
        r_code += '''#MS2 id is %s\n''' % fragment

        rt_list = ','.join(map(str, display_data[tg][sample]['ms2']['rt_list'][fragment]))
        r_code += '''rt = c(%s)\n''' % rt_list
        i_list = ','.join(map(str, [round(float(x) * 100.0 / max_intensity_ms2, 2) for x in display_data[tg][sample]['ms2']['i_list'][fragment]]))
        r_code += '''int = c(%s)\n''' % i_list
        r_code += '''lines(rt, int, type = "l", col = %d, lwd = %.1f)\n''' % (
            transition_color_code, parameters.PLOT_LINE_WIDTH)
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

