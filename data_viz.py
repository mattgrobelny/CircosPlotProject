
###############################################################################
###############################################################################

# create a dictionary for chromosome size
chr_size_dic = {'groupI': 28185914, 'groupII': 23295652, 'groupIII': 16798506, 'groupIV': 32632948,
                'groupIX': 20249479, 'groupV': 12251397, 'groupVI': 17083675, 'groupVII': 27937443,
                'groupVIII': 19368704, 'groupX': 15657440, 'groupXI': 16706052, 'groupXII': 18401067,
                'groupXIII': 20083130, 'groupXIV': 15246461, 'groupXIX': 20240660, 'groupXV': 16198764,
                'groupXVI': 18115788, 'groupXVII': 14603141, 'groupXVIII': 16282716, 'groupXX': 19732071,
                'groupXXI': 11717487}

chrm_name_order_list = ['groupI', 'groupII', 'groupIII', 'groupIV',
                        'groupV', 'groupVI', 'groupVII',
                        'groupVIII', 'groupIX', 'groupX', 'groupXI', 'groupXII',
                        'groupXIII', 'groupXIV', 'groupXV',
                        'groupXVI', 'groupXVII', 'groupXVIII', 'groupXIX', 'groupXX', 'groupXXI']

stats_dic = {'Fst': (0, 1), 'Div': (0, 1)}

stat_list = ['Fst', 'Div']

color_grad_dic = {'Fst': 'Greens', 'Div': 'seismic'}


###############################################################################
viz_parameters = {'total_genome_size': sum(chr_size_dic.values()),
                  'number_of_chr': len(chr_size_dic),
                  'rad_inner': 250,
                  'ring_gap': 10,
                  'arc_padding_in_degrees': 2,
                  'last_degree_end': 0,
                  'ring_width': 35,
                  'total_degrees': 0,
                  '10mb_step_off_set': 32,
                  'font_size': 5,  # need to test sizes
                  'width': 2,   # must be float ? need to test
                  # a sequence specifying alternate lengths of on and off
                  # stroke portions.
                  'dash_pattern': [3, 1],
                  'label_units': "Mb",
                  'key_loc_offset': 90,
                  'key_width': 30,
                  'key_height': 120,
                  'key_sep_distance': 14,
                  'key_degree_off_set': -10,
                  "key_label_font_x": 0.7
                  }

# calculate number of degrees per nucleotide
viz_parameters['degree_per_nuc'] = float(360 - (viz_parameters['number_of_chr'] * viz_parameters[
                                         'arc_padding_in_degrees'])) / float(viz_parameters['total_genome_size'])

img = {}
img['height'] = 800
img['width'] = 800
img['center_x'] = img['width'] / 2.0
img['center_y'] = img['height'] / 2.0
img['font_size'] = 16

###############################################################################
# available colors specturms names
cmaps = [('Perceptually Uniform Sequential',
          ['viridis', 'inferno', 'plasma', 'magma']),
         ('Sequential',     ['Blues', 'BuGn', 'BuPu',
                             'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                             'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                             'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
         ('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool',
                             'copper', 'gist_heat', 'gray', 'hot',
                             'pink', 'spring', 'summer', 'winter']),
         ('Diverging',      ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                             'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                             'seismic']),
         ('Qualitative',    ['Accent', 'Dark2', 'Paired', 'Pastel1',
                             'Pastel2', 'Set1', 'Set2', 'Set3']),
         ('Miscellaneous',  ['gist_earth', 'terrain', 'ocean', 'gist_stern',
                             'brg', 'CMRmap', 'cubehelix',
                             'gnuplot', 'gnuplot2', 'gist_ncar',
                             'nipy_spectral', 'jet', 'rainbow',
                             'gist_rainbow', 'hsv', 'flag', 'prism'])]
# http://matplotlib.org/examples/color/colormaps_reference.html

color_stat_mapper_dic = {}

for stat in stat_list:
    color_stat_mapper_dic[stat] = {}

import cairocffi as cairo
# import cairo
import sys
import getopt
import math
import numpy as np

import matplotlib.pyplot as plt
import circos_draw

###############################################################################
#
# Define the path our file
#
#outpath = '/home/a-m/ib501_stud12/shell/data_viz/data_viz.pdf'
#outpath = '/home/mgrobelny/Scripts/github/Data-viz-Circle-plot/data_viz.pdf'
outpath = '/Users/matt/github/Data-viz-Circle-plot/data_viz.pdf'

ps = cairo.PDFSurface(str(outpath), float(img['height']), float(img['width']))
cr = cairo.Context(ps)

# import drawing functions

###############################################################################
# # default parameters

# argv = sys.argv[1:]
# try:
#     opts, args = getopt.getopt(argv, "hs:")
# except getopt.GetoptError:
#     print 'kspec.py -k <kmer_size> -x <x_axis_max> -t <type>[fasta|fastq] -f <inputfile>'
#     sys.exit(2)
# for opt, arg in opts:
#     if opt == '-h':
#         print "#--- K-mer frequency graphing script ---#\n"
#         print "Usage:"
#         print 'kspec.py -k <kmer_size> -x <x_axis_max> -t <type>[fasta|fastq] -f <inputfile> \n'
#         print "Goals:"
#         print "1) Take in fastq file and kmerize it and output kmer occurence frequnecy"
#         print "2) Output graph of kmer occurence frequnecy"
#         print "3) Output kmer occurence frequnecy to .tsv file"
#         print "\n"
#
#         sys.exit()
#     elif opt in ("-k"):
#         kmer = arg
#     elif opt in ("-x"):
#         xmax = arg
#     elif opt in ("-f"):
#         file_name = arg
#     elif opt in ("-t"):
#         file_type = arg
# print "Input file:", file_name
# print "Input file type:", file_type
# print "Kmer size:", kmer
# print "X-axis max kmer count:", xmax
# print " "
#
#
# ###############################################################################


###############################################################################
# # Data import

# Create fst data dictionary
fst_stats = {}
rna_stats = {}
for chrm_name in chrm_name_order_list:
    fst_stats[chrm_name] = []
    rna_stats[chrm_name] = []

level_to_dic = {0: fst_stats,
                1: rna_stats}

# Add each chromosome to the dictionary and store the
# basepair and statistical value
# add each chrm, bp and stat pt to dictionary

# Open fst Data file
in_file = '/Users/matt/github/Data-viz-Circle-plot/Pop_fst_out.tsv'
fh1_fst_file = open(in_file, 'r')

# skip header
next(fh1_fst_file)
for line in fh1_fst_file:
    # strip new line char
    line = line.strip('\n')
    # remove spaces
    line = line.replace(" ", "")
    # split tabs
    line = line.split('\t')
    # append data to each dictionary of  list
    fst_stats[line[0]].append([line[1], line[2]])
fh1_fst_file.close

# repeate for Div data_viz
in_file2 = '/Users/matt/github/Data-viz-Circle-plot/Pop_div_data.tsv'
fh2_Div_file = open(in_file2, 'r')

# skip header
next(fh2_Div_file)
for line in fh2_Div_file:
    # strip new line char
    line = line.strip('\n')
    # remove spaces
    line = line.replace(" ", "")
    # split tabs
    line = line.split('\t')
    # append data to each dictionary of list
    rna_stats[line[0]].append([line[1], line[2]])
fh2_Div_file.close

#  Thus:
#  For rna_stats['chrmII'] outputs ["bp","stat"]
#  rna_stats['chrmII'][0] = ["bp","stat"]
#  rna_stats['chrmII'][0][0] = Base pair
#  rna_stats['chrmII'][0][1] = stats

###############################################################################
# Data normalization
for chrm_name in chrm_name_order_list:
    #data_norm(chrm_name,fst_stats,fst_stats[chrm_name], "norm")
    circos_draw.data_norm(chrm_name, rna_stats, rna_stats[chrm_name], "log")


###############################################################################
# Map colors to stat

circos_draw.stat_to_color('Fst', "norm", 'False')
circos_draw.stat_to_color('Div', "norm", "False")


###############################################################################
# Draw final image
circos_draw.draw_chrom_arc_w_label(chrm_name_order_list, 2, 1, 1, "def")
###############################################################################

#
# End of file
#
# Close the file
#
cr.show_page()
