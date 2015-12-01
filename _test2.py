__author__ = 'Tiannan Guo, ETH Zurich 2015'

import re
with open('_test2.txt', 'w') as myfile, open('debug_png_id_200.R', 'r') as infile:
    for line in infile:
        m = re.search('plot', line)
        if m:
            myfile.write(line)


