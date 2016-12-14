#!/usr/bin/python

#import cairocffi as cairo
import cairo
import sys
import getopt
import math


###############################################################################
###############################################################################

# create a dictionary for chromosome size
chr_size_dic = {'groupI' : 28185914, 'groupII'   : 23295652, 'groupIII'   : 16798506, 'groupIV'  : 32632948,
'groupIX'   : 20249479, 'groupV'    : 12251397, 'groupVI'    : 17083675, 'groupVII' : 27937443,
'groupVIII' : 19368704, 'groupX'    : 15657440, 'groupXI'    : 16706052, 'groupXII' : 18401067,
'groupXIII' : 20083130, 'groupXIV'  : 15246461, 'groupXIX'   : 20240660, 'groupXV'  : 16198764,
'groupXVI'  : 18115788, 'groupXVII' : 14603141, 'groupXVIII' : 16282716, 'groupXX'  : 19732071,
'groupXXI'  : 11717487}

chrm_name_order_list= ('groupI', 'groupII','groupIII', 'groupIV',
'groupIX', 'groupV', 'groupVI', 'groupVII',
'groupVIII', 'groupX', 'groupXI', 'groupXII',
'groupXIII', 'groupXIV', 'groupXIX', 'groupXV',
'groupXVI', 'groupXVII', 'groupXVIII', 'groupXX',
'groupXXI')
###############################################################################
viz_parameters = {'total_genome_size': sum(chr_size_dic.values()),
'number_of_chr': len(chr_size_dic),
'rad_inner' : 250,
'ring_gap': 10,
'arc_padding_in_degrees': 2,
'last_degree_end': 0,
'ring_width': 35
'total_degrees': 0,
'10mb_step_off_set': 5,
'width': 0.5   # must be float ? need to test
#'fill_color' :'0.4,0.4,0.4' ,
#'trim_color' : '0,0,0'
}

viz_parameters['degree_per_nuc'] = float(360 - (viz_parameters['number_of_chr'] * viz_parameters['arc_padding_in_degrees'])) / float(viz_parameters['total_genome_size'])

img = {}
img['height']     = 800
img['width']      = 800
img['center_x']   = img['width'] / 2.0
img['center_y']   = img['height'] / 2.0
img['font_size']  = 16

###############################################################################
#
# Define the path our file
#
#outpath = '/home/a-m/ib501_stud12/shell/data_viz/data_viz.pdf'
outpath='/home/mgrobelny/Scripts/github/Data-viz-Circle-plot/data_viz.pdf'
ps = cairo.PDFSurface(str(outpath), float(img['height']), float(img['width']))
cr = cairo.Context(ps)


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
# # Progress bar is not my own work from:
# # https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
# #
# def progress(count, total, suffix=''):
#     bar_len = 60
#     filled_len = int(round(bar_len * count / float(total)))
#
#     percents = round(100.0 * count / float(total), 1)
#     bar = '=' * filled_len + '-' * (bar_len - filled_len)
#
#     sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
#     sys.stdout.flush()
# ##################################################

###############################################################################
#
# Convert a radius and a span of degrees into X, Y coordinates #
def get_x_y_coordinates(center_x, center_y, degree, radius):
    if degree <= 90:
        theta = float(degree)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x + adj_side
        y = center_y + opp_side
    elif degree <= 180:
        theta = float(degree - 90.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x - opp_side
        y = center_y + adj_side
    elif degree <= 270:
        theta = float(degree - 180.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x - adj_side
        y = center_y - opp_side

    else:
        theta = float(degree - 270.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x + opp_side
        y = center_y - adj_side
    return (x, y)

# int to roman not my code from:
# https://www.safaribooksonline.com/library/view/python-cookbook/0596001673/ch03s24.html
def int_to_roman(input):
    """ Convert an integer to a Roman numeral. """
    if not isinstance(input, type(1)):
        raise TypeError, "expected integer, got %s" % type(input)
    if not 0 < input < 4000:
        raise ValueError, "Argument must be between 1 and 3999"
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
    result = []
    for i in range(len(ints)):
        count = int(input / ints[i])
        result.append(nums[i] * count)
        input -= ints[i] * count
    return ''.join(result)

###############################################################################

# Arc drawing functions

# Draw a circle of arcs based on a list of chr which corresponed to the chrm size dic
def draw_label(text, x, y):

    # Font
    cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    # Set the font size
    cr.set_font_size(font_size)
    # font color
    cr.set_source_rgb(0, 0, 0)

    # Get the size of the text we want to write, returns a tuple:
    #   (x, y, width, height, dx, dy)
    #
    textents = cr.text_extents(text)
    text_width = textents[2]
    text_height = textents[3]
    #
    # Where you want to draw text may need to be adjusted,
    # depending on the size of the text.
    cr.move_to(x, y)
    cr.show_text(text)

def chrm_arc(chrm_name, level, trim):

    # Create intial arc
    radius = viz_parameters['rad_inner'] + viz_parameters['ring_gap'] * level + viz_parameters['ring_width'] * level

    # calculate the number of degrees the are will span based on chrm size
    viz_parameters['total_degrees'] = float(viz_parameters['degree_per_nuc']) * float(chr_size_dic[chrm_name])

    # draw first arc based on the ending of the pervious arc
    sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], viz_parameters['last_degree_end'], radius)

    cr.move_to(sx, sy)
    start_deg = viz_parameters['last_degree_end']
    end_deg = viz_parameters['last_degree_end'] + viz_parameters['total_degrees']

    # draw outer arc
    cr.new_sub_path()
    cr.arc_negative(img['center_x'], img['center_y'], radius, math.radians(end_deg),math.radians(start_deg))

    # draw line to inner arc
    sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], viz_parameters['last_degree_end'], radius - viz_parameters['ring_width'])
    cr.line_to(sx, sy)

    # draw reverse arc
    cr.arc(img['center_x'], img['center_y'], radius - viz_parameters['ring_width'], math.radians(start_deg), math.radians(end_deg))

    cr.close_path()

    ############################
    # Add black trim to chrm arcs 0 no , 1 yes
    if trim == 0:
        # cr.set_source_rgb(viz_parameters['fill_color'])

        # fill with grey color
        cr.set_source_rgb(0.4, 0.4, 0.4)
        cr.fill()
    else:
        #cr.set_source_rgb(viz_parameters['trim_color'])

        #  trim in black
        cr.set_source_rgb(0, 0, 0)
        cr.stroke()

    ############################

    # Update the end of viz parameter[last_degree_end] + padding --> for next arc start degree
    viz_parameters['last_degree_end'] = float(viz_parameters['last_degree_end']) + float(viz_parameters['total_degrees']) + float(viz_parameters['arc_padding_in_degrees'])

def draw_10mb_labels(chrm_name, level):
    # Add 10mb labels
    ten_mb_step = 0
    ten_mb = int(chr_size_dic[chrm_name] / 1000000)
    for i in range(ten_mb):
        # caculate 10mb step
        ten_mb_step = ten_mb_step + i * 1000000

        # determine degree of 10mb step line
        ten_mb_step_degree = float(ten_mb_step * viz_parameters['degree_per_nuc'] + viz_parameters['last_degree_end'])

        # find the x and y pos of the location of the 10mb step
        sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], ten_mb_step_degree, radius - viz_parameters['10mb_step_off_set'])

        # move to that location
        cr.move_to(sx, sy)

        # find the end of the line
        end_of_line = radius + viz_parameters['ring_gap'] * (level) + viz_parameters['ring_width'] * level + viz_parameters['10mb_step_off_set']
        sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], ten_mb_step_degree, end_of_line)
        cr.line_to(sx, sy)

        if i % 2 == 0:
            # darker grey line color
            cr.set_source_rgb(0.3, 0.3, 0.3)

            # stroke a thicker line
            set_line_width(viz_parameters['width'] + 0.1)

        else:
            # ligher grey line color
            cr.set_source_rgb(0.4, 0.4, 0.4)

            # stroke a thinner line
            set_line_width(viz_parameters['width'])

            set_line_width(width)
            dash = 0

        cr.stroke


def draw_chrom_arc(chrm_list, level, trim):
    for key in chrm_list:
        chrm_arc(key, level, trim)
    viz_parameters['last_degree_end'] = 0


# Test1
for i in range(3):
    draw_chrom_arc(chrm_name_order_list, i, 0)
    draw_chrom_arc(chrm_name_order_list, i, 1)


###############################################################################

# # Data import
# in_file =
# fh1 = open(in_file, 'r')
# fst_stats = {}
# # Add each chromosome to the dictionary and store the
# # basepair and statistical value
#
# fst_stats[chr]['bp']   = []
# fst_stats[chr]['stat'] = []
# rna_stats = {}
# rna_stats[chr]['bp']   = []
# ran_stats[chr]['stat'] = []

###############################################################################
# End of file
#
# Close the file
#
cr.show_page()
