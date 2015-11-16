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

