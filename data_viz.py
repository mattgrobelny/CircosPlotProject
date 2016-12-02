# Cairo drawing module
#
import cairocffi as cairo

###############################################################################
#
# The img dictionary will hold image-specific dimensions.
#
img = {}
img['height']     = 800
img['width']      = 800
img['center_x']   = img['width']  / 2.0
img['center_y']   = img['height'] / 2.0
img['font_size']  = 16;

###############################################################################
#
# Define the path our file
#
outpath = "/home/a-m/ib501_stud12/shell/Final_hw"
ps = cairo.PDFSurface(outpath, img['height'], img['width'])
cr = cairo.Context(ps)


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
        y = centre_y - opp_side

    else:
        theta = float(degree - 270.0)
        opp_side = radius * math.sin(math.radians(theta))
        adj_side = radius * math.cos(math.radians(theta))
        x = center_x + opp_side
        y = center_y - adj_side
    return (x, y)

###############################################################################
#
# Choose a font
#
cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
# Set the font size
cr.set_font_size(font_size)
# Choose a font color
cr.set_source_rgb(red, green, blue)
#
# Get the size of the text we want to write, returns a tuple:
#   (x, y, width, height, dx, dy)
#
textents = cr.text_extents(text)
text_width = textents[2]
text_height = textents[3]
#
# Where you want to draw text may need to be adjusted,
# depending on the size of the text.
#
# cr.move_to(x, y)
#
# Finally draw the text.
#
# cr.show_text(text)

###############################################################################

# create a dictionary for chromosome size
chr_size_dic = {'groupI'    : 28185914, 'groupII'   : 23295652, 'groupIII'   : 16798506, 'groupIV'  : 32632948,
'groupIX'   : 20249479, 'groupV'    : 12251397, 'groupVI'    : 17083675, 'groupVII' : 27937443,
'groupVIII' : 19368704, 'groupX'    : 15657440, 'groupXI'    : 16706052, 'groupXII' : 18401067,
'groupXIII' : 20083130, 'groupXIV'  : 15246461, 'groupXIX'   : 20240660, 'groupXV'  : 16198764,
'groupXVI'  : 18115788, 'groupXVII' : 14603141, 'groupXVIII' : 16282716, 'groupXX'  : 19732071,
'groupXXI'  : 11717487}

###############################################################################
viz_parameters = {'total_genome_size': int(sum(chr_size_dic.items())),
'number_of_chr': 21,
'degree_per_nuc': float((360 - int(viz_parameters['number_of_chr']))/ viz_parameters['total_genome_size']),
'rad_inner' : 250,
'rad_outer': 300,
'ring_gap': 10,
'arc_padding_in_degrees': 1,
'last_degree_end': 0,
'ring_width': 25}

def chrm_arc(chrm_name, rad_type):
    if rad_type == "inner":
        radius = viz_parameters['rad_inner']
    else:
        radius = viz_parameters['rad_outer']

    total_degrees = (viz_parameters['degree_per_nuc']) * chr_size_dic[chrm_name]
    sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], viz_parameters['last_degree_end'], radius)
    cr.move_to(sx, sy)
    start_deg = viz_parameters['last_degree_end']
    end_deg = viz_parameters['last_degree_end'] + total_degrees

    # draw outer arc
    cr.set_source_rgb(0, 0, 0)
    cr.arc(img['center_x'], img['center_y'], radius, math.radians(start_deg), math.radians(end_deg))

    # draw line to inner arc
    sx, sy = get_x_y_coordinates(img['center_x'], img['center_y'], viz_parameters['last_degree_end'], radius - viz_parameters['ring_width'])
    cr.line_to(sx, sy)

    # draw reverse arc
    cr.arc_negative(img['center_x'], img['center_y'], radius - viz_parameters['ring_width'], math.radians(start_deg), math.radians(end_deg))

    # close arc
    cr.close_path()

    # write  line
    cr.stroke()

    # fill with gray
    cr.set_source_rgb(0.4, 0, 0)
    cr.fill()

    # Update the end of viz parameter[last_degree_end] + padding --> for next arc start degree
    viz_parameters['last_degree_end'] = viz_parameters['last_degree_end'] + total_degrees + viz_parameters['arc_padding_in_degrees']

# test1
chrm_arc('groupI', 'inner')

###############################################################################
###############################################################################
# End of file
#
# Close the file
#
cr.show_page()
