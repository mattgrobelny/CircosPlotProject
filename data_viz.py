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
'groupV', 'groupVI', 'groupVII',
'groupVIII', 'groupIX', 'groupX', 'groupXI', 'groupXII',
'groupXIII', 'groupXIV', 'groupXV',
'groupXVI', 'groupXVII', 'groupXVIII', 'groupXIX','groupXX','groupXXI')

stats_dic= {'Div' :(0,1),
'Fst': (0, 1),
'Rand': (0,1)
}
stat_list = ('Div','Fst', 'Rand')
color_grad_dic= {'Div':[(0.2,0.2,0.6),(0.5,0.6,0.7)],
'Fst':  [(0.3,0.2,0.6),(0.5,0.6,0.7)],
'Rand': [(0.2,0.4,0.6),(0.5,0.9,0.7)]

}
###############################################################################
viz_parameters = {'total_genome_size': sum(chr_size_dic.values()),
'number_of_chr': len(chr_size_dic),
'rad_inner' : 250,
'ring_gap': 10,
'arc_padding_in_degrees': 2,
'last_degree_end': 0,
'ring_width': 35,
'total_degrees': 0,
'10mb_step_off_set':32,
'font_size': 5, # need to test sizes
'width': 2,   # must be float ? need to test
'dash_pattern': [3,1], # a sequence specifying alternate lengths of on and off stroke portions.
'label_units': "Mb",
'key_loc_offset' : 90,
'key_width': 20,
'key_height' : 100,
'key_sep_distance':14,
'key_degree_off_set': -10,
"key_label_font_x": 0.5

#'fill_color' :'0.4,0.4,0.4' ,
#'trim_color' : '0,0,0'
}

# calculate number of degrees per nucleotide
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
#

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

# --------  Drawing functions -------- #

# Draw a circle of arcs based on a list of chr which corresponed to the chrm size dic
def draw_label(text, x, y, font_size,working_degree):

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
    if working_degree >=100 and working_degree <= 280:
        centered_x= x
        centered_y = y + text_height/2
    else:
        centered_x= x - text_width
        centered_y =y + 1

    cr.move_to(centered_x, centered_y)
    # cr.move_to(x,y)
    cr.show_text(text)

# Draw 10mb label markers
def draw_10mb_labels(chrm_list, level):

    for chrm_name_it in chrm_list:

        break_size = 5000000 #bases
        # find how many 10mb breaks there are for chrm_name
        five_mb_break = int(chr_size_dic[chrm_name_it] / break_size)

        # determine degree of 10mb step line
        five_mb_step_degree = float(break_size * viz_parameters['degree_per_nuc'])

        # for i number of breaks draw a line every 5mb and a
        for i in range(1,five_mb_break+1):

            # calculate 5mb step
            working_degree = five_mb_step_degree * i + viz_parameters['last_degree_end']

            # find the x and y pos of the location of the 10mb step
            sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], working_degree, viz_parameters['rad_inner']  - viz_parameters['10mb_step_off_set']* 1.4)

            # move to that location
            cr.move_to(sx, sy)

            # find the end of the line
            end_of_line = viz_parameters['rad_inner'] + viz_parameters['ring_gap'] * level + viz_parameters['ring_width'] * (level -1)
            sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], working_degree, end_of_line)
            # write line
            cr.line_to(sx, sy)

            # made dashed line
            cr.set_dash(viz_parameters['dash_pattern'])

            if i % 2 == 0:
                # darker grey line color
                cr.set_source_rgb(0.2, 0.2, 0.2)

                # stroke a thicker line
                cr.set_line_width(viz_parameters['width'] + 0.1)

                cr.stroke()

                # find the x and y location for the 10m label
                label_x, label_y = get_x_y_coordinates(img['center_x'], img['center_y'], working_degree, viz_parameters['rad_inner']  - viz_parameters['10mb_step_off_set']* 1.49)
                #working_degree, viz_parameters['rad_inner'] - viz_parameters['10mb_step_off_set']*1.75)

                # write label name w/ units
                label = str(int(i * 5)) + viz_parameters['label_units']

                # pass to draw label function
                draw_label(label, label_x, label_y, 7, working_degree)

            else:
                # ligher grey line color
                cr.set_source_rgb(0.4, 0.4, 0.4)

                # stroke a thinner line
                cr.set_line_width(viz_parameters['width']-1)
                cr.stroke()

        # Update where the start of next chrm is os labeling can be indexed corretly
        viz_parameters['total_degrees'] = float(viz_parameters['degree_per_nuc']) * float(chr_size_dic[chrm_name_it])
        viz_parameters['last_degree_end'] = float(viz_parameters['last_degree_end']) + float(viz_parameters['total_degrees']) + float(viz_parameters['arc_padding_in_degrees'])

# Draw chrm name labels use ither provided labels from list or generate new names with roman number (roman= 1)
def chrm_label(chrm_list, total_levels, roman):
    count = 0
    for chrm_name in chrm_list:
        degree_for_label = float(viz_parameters['degree_per_nuc']) * float(chr_size_dic[chrm_name])/2
        working_degree = viz_parameters['last_degree_end'] + degree_for_label
        radian_for_label = viz_parameters['rad_inner'] + viz_parameters['ring_gap'] * total_levels + viz_parameters['ring_width'] * total_levels

        label_x, label_y = get_x_y_coordinates(img['center_x'], img['center_y'], working_degree, radian_for_label)
        #working_degree, viz_parameters['rad_inner'] - viz_parameters['10mb_step_off_set']*1.75)
        label = ""

        if roman == 1:
            label = str(int_to_roman(count + 1))
            count = count +1
        else:
            label = str(chrm_name)

        # pass to draw label function
        draw_label(label, label_x, label_y, 12, working_degree)


        # Update where the start of next chrm is os labeling can be indexed corretly
        viz_parameters['total_degrees'] = float(viz_parameters['degree_per_nuc']) * float(chr_size_dic[chrm_name])
        viz_parameters['last_degree_end'] = float(viz_parameters['last_degree_end']) + float(viz_parameters['total_degrees']) + float(viz_parameters['arc_padding_in_degrees'])

def color_key(total_levels,location,trim): #min, max,color_start, color_end,
    if location == "top_left":
        working_degree = 225
    elif location == "bottom_right":
        working_degree = 45
    elif location == "bottom_left":
        working_degree = 135
    else:
        # default to top_right location for key
        working_degree = 315

    radius = viz_parameters['rad_inner'] + viz_parameters['ring_gap'] * total_levels + viz_parameters['ring_width'] * total_levels + viz_parameters['key_loc_offset']
    sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], working_degree + viz_parameters['key_degree_off_set'], radius)
    sx_key = 0
    for i in range(total_levels):
        sx_key = sx + (viz_parameters['key_width'])* i +viz_parameters['key_sep_distance'] * i
        cr.move_to(sx_key, sy)
        cr.rectangle(sx_key, sy, viz_parameters['key_width'], viz_parameters['key_height'])
        cr.close_path()
        cr.set_line_width(viz_parameters['width'] - 1)
        cr.set_dash([])

        ############################
        # Add black trim to key 0 no , 1 yes

        if trim == 0:
            # fill with color gradiant
            x0 = sx_key +viz_parameters['key_width']/2
            y0 = sy
            x1 = x0
            y1 = sy - viz_parameters['key_height']
            # grad_fil = cr.LinearGradient(x0, y0, x1,y1)
            # grad_fil.add_color_stop_rgba(0, "%s"%(color_grad_dic[stat_list[i][0]]), 1)
            # grad_fil.add_color_stop_rgba(0, "%s"%(color_grad_dic[stat_list[i][0]]), 1)

            lg1 = cr.LinearGradient(x0, y0, x1,y1)

            count = 1

            a = 0.1
            while a < 1.0:
                if count % 2:
                    lg1.add_color_stop_rgba(a, 0, 0, 0, 1)
                else:
                    lg1.add_color_stop_rgba(a, 1, 0, 0, 1)
                a = a + 0.1
                count = count + 1
            cr.set_source_rgb(lg1)
            cr.fill()
            cr.stroke()
        else:
            #cr.set_source_rgb(viz_parameters['trim_color'])

            #  trim in black
            cr.set_source_rgb(0, 0, 0)
            cr.stroke()

        # Draw Key labels

        # Draw stat name above key
        draw_label(stat_list[i], sx_key + viz_parameters['key_width']/2 +4, sy-3, 12* viz_parameters['key_label_font_x'],working_degree)

        # Draw min labels
        draw_label(str(stats_dic[stat_list[i]][0]), sx_key-2, sy + viz_parameters['key_height'], 9* viz_parameters['key_label_font_x'],working_degree)

        # Draw max label
        draw_label(str(stats_dic[stat_list[i]][1]), sx_key- 2, sy, 9* viz_parameters['key_label_font_x'],working_degree)


# Draw chrm arc for a given level w (1) or wo (0) balck trim
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
    cr.set_line_width(viz_parameters['width'] - 1)
    cr.set_dash([])
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

# Draw all chrm arc for a given level w (1) or wo (0) balck trim
def draw_chrom_arc(chrm_list, level, trim):
    for key in chrm_list:
        chrm_arc(key, level, trim)
    viz_parameters['last_degree_end'] = 0

# -------- Main Drawing function -------- #

# Draw all chrm arc for a given level w (1) or wo (0) balck trim
# w 10mb labels
def draw_chrom_arc_w_label(chrm_list, total_levels, trim, roman, location):
    viz_parameters['last_degree_end'] = 0

    # draw all breaks and labels
    draw_10mb_labels(chrm_list, total_levels)
    viz_parameters['last_degree_end'] = 0
    # draw chm labels
    chrm_label(chrm_list, total_levels, roman)
    viz_parameters['last_degree_end'] = 0

    #Draw all chrm arcs
    if trim == 0:
        for i in range(total_levels):
            draw_chrom_arc(chrm_list, i, 0)
    else:
        for i in range(total_levels):
            draw_chrom_arc(chrm_list, i, 0)
            draw_chrom_arc(chrm_list, i, 1)

    # Draw Key
    color_key(total_levels, location, trim)
###############################################################################
# Test 2 - should output

draw_chrom_arc_w_label(chrm_name_order_list, 3, 1, 1,"def")
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
