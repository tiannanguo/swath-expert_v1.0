__author__ = 'Tiannan Guo, ETH Zurich 2015'

# parameters
PEAK_TOLERANCE = 6  # seconds, use to find if the peak is found
title = ['transition_name','transition_group_id', 'best_rt', 'best_sample', 'best_score',
         'protein', 'irt', 'precursor_mz', 'transition_mz', 'transition_i']
MIN_FRAGMENTS = 4
MAX_RT_TOLERANCE = 10
MAX_PEAK_WIDTH = 100
PEAK_WIDTH_FOLD_VARIATION = 2.0
PEAK_SHAPE_FOLD_VARIATION = 3.0
#MIN_TRANSITION_NUMBER = 5  # after optimization, I find set min 5 transitions is good

PLOT_LINE_WIDTH = 2.0  # the thickness of chrom curve

# for golden 30 samples
# figure_num_per_row = 10
# figures_num_rows = 3
# PNG_FILE_WIDTH = figure_num_per_row * 15
# PNG_FILE_HEIGHT = figures_num_rows * 600

# for golden 90 samples
figure_num_per_row = 30
figures_num_rows = 3
PNG_FILE_WIDTH = figure_num_per_row * 150
PNG_FILE_HEIGHT = figures_num_rows * 600



title_font_size = 2
plot_out_margin_area_south = 10
plot_out_margin_area_west = 5
plot_out_margin_area_north = 15
plot_out_margin_area_east = 5
subplot_out_margin_area_south = 5
subplot_out_margin_area_west = 2
subplot_out_margin_area_north = 5
subplot_out_margin_area_east = 2
frame_of_plot = 'TRUE' # whether show a frame outside the subplots

