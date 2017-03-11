#!/usr/bin/python

import sys

PRECISION = 1

def help():
    print("usage: {0} [input_file] [output_file]".format(sys.argv[0]))
    print("input_file format:")
    print("------------------")
    print("[PRECISION]")
    print("[X SPACES] [Y SPACES]")
    print("[minX RANGE] [maxX RANGE]")
    print("[minY RANGE] [maxY RANGE]")
    print("[X Y DATA POINT SEPEARTED BY A SPACE FOR EVERY LINE]")

def read_file(fname):
    global PRECISION
    with open(fname) as f:
        PRECISION = [int(x) for x in next(f).split()][0]
        xspaces, yspaces = [int(x) for x in next(f).split()]
        xmin, xmax = [float(x) for x in next(f).split()]
        ymin, ymax = [float(x) for x in next(f).split()]
        xdata = []
        ydata = []
        for line in f:
            x, y = [float(a) for a in line.split()]
            xdata.append(x)
            ydata.append(y)
    return [xspaces, yspaces, [xmin, xmax], [ymin, ymax], xdata, ydata]

def write_file(fname, scales, coords, datapoints, gradient, intersect, fitcoords, gravitypoint):
    with open(fname, 'w+') as f:
        f.write("SCALES:\n")
        f.write("X: {0}\t Xmm: {1}\n".format(scales['x']['scale'], scales['x']['scale_mm']))
        f.write("Y: {0}\t Ymm: {1}\n".format(scales['y']['scale'], scales['y']['scale_mm']))

        f.write("\nCOORDINATES:\n")
        f.write("X: {0}\n".format(", ".join(map(str, coords['x']))))
        f.write("Y: {0}\n".format(", ".join(map(str, coords['y']))))

        f.write("\nDATAPOINTS:\n")
        for i in datapoints:
            f.write("X: {0}\nXoffset: {1}\nY: {2}\nYoffset: {3}\n\n".format(i['x'][0], i['x'][1], i['y'][0], i['y'][1]))

        f.write("FITLINE:\n")
        f.write("Gradient: {0}\nIntersect: {1}\n".format(gradient, intersect))
        f.write("\nFITPOINTS:\n")
        for i in fitcoords:
            f.write("X: {0}\nY: {1}\nYoffset: {2}\n\n".format(i['x'][0], i['y'][0], i['y'][1]))

        f.write("GRAVITYPOINT:\n")
        f.write("X: {0}\nXoffset: {1}\nY: {2}\nYoffset: {3}\n\n".format(gravitypoint['x'][0], gravitypoint['x'][1], gravitypoint['y'][0], gravitypoint['y'][1]))
        f.write("X: {0}\nY: {1}\n".format(gravitypoint['x'][2], gravitypoint['y'][2]))

def find_mean(vals):
    total = 0.0
    for i in vals:
        total += i
    return total / len(vals)

def find_gradient(xvals, yvals, xmean, ymean):
    num = 0.0
    den = 0.0
    for i in range(0, len(xvals)):
        xdiff = xvals[i] - xmean
        num += (xdiff) * (yvals[i] - ymean)
        den += xdiff ** 2
    return num / den

def find_scale(rang, spaces):
    vmax = rang[1]
    vmin = rang[0]

    sv = round((vmax - vmin) / spaces, PRECISION)
    sv_small = round(sv / 10, PRECISION + 1)

    return {'scale': sv, 'scale_mm': sv_small, 'min': vmin, 'spaces': spaces}

def find_coords(scales):
    arr = []
    for i in range(scales['spaces']):
        arr.append(round(scales['scale'] * i + scales['min'], PRECISION))
    return arr

def find_nearest_coord(xval, yval, xscales, yscales, xcoords, ycoords):
    # calculate x coordinate
    xcoord = xcoords[len(xcoords) - 1]
    for i in reversed(xcoords):
        if xval - i >= 0:
            xcoord = i
            break

    # calculate y coordinate
    ycoord = ycoords[len(ycoords) - 1]
    for i in reversed(ycoords):
        if yval - i >= 0:
            ycoord = i
            break

    # calculate offsets
    xoffset = xval - xcoord
    yoffset = yval - ycoord

    xoffseti = round(xoffset / xscales['scale_mm'])
    yoffseti = round(yoffset / yscales['scale_mm'])

    if xoffseti > 10:

        xcoord = round(xcoord + xscales['scale'], PRECISION)
        xoffseti = 0

    if yoffseti > 10:
        ycoord = round(ycoord + yscales['scale'], PRECISION)
        yoffseti = 0

    return {'x': [xcoord, xoffseti], 'y': [ycoord, yoffseti]}

def generate_fitline_coords(xdata, ydata, xscales, yscales, xcoords, ycoords):

    xmean = find_mean(xdata)
    ymean = find_mean(ydata)

    m = find_gradient(xdata, ydata, xmean, ymean)

    #m = ymean / xmean
    b = ymean - m * xmean

    coords = []
    
    yval = round(m * xcoords[0] + b, PRECISION)
    coords.append(find_nearest_coord(xcoords[0], yval, xscales, yscales, xcoords, ycoords))

    yval = round(m * xcoords[len(xcoords) - 1] + b, PRECISION)
    coords.append(find_nearest_coord(xcoords[len(xcoords) - 1], yval, xscales, yscales, xcoords, ycoords))

    return [m, b, coords]

def find_gravity_point(xdata, ydata, xscales, yscales, xcoords, ycoords):
    xmean = find_mean(xdata)
    ymean = find_mean(ydata)

    r = find_nearest_coord(xmean, ymean, xscales, yscales, xcoords, ycoords)
    r['x'].append(round(xmean, PRECISION))
    r['y'].append(round(ymean, PRECISION))

    return r

if __name__ == "__main__":

    if len(sys.argv) < 3:
        help()
        sys.exit(0)

    xspaces, yspaces, xrang, yrang, xdata, ydata = read_file(sys.argv[1])

    xscales = find_scale(xrang, xspaces)
    yscales = find_scale(yrang, yspaces)

    xcoords = find_coords(xscales)
    ycoords = find_coords(yscales)

    datapoints = []
    for i in range(0, len(xdata)):
        datapoints.append(find_nearest_coord(xdata[i], ydata[i], xscales, yscales, xcoords, ycoords))

    gradient, intersect, fitcoords = generate_fitline_coords(xdata, ydata, xscales, yscales, xcoords, ycoords)

    gravitypoint = find_gravity_point(xdata, ydata, xscales, yscales, xcoords, ycoords)

    scales = {'x': xscales, 'y': yscales}
    coords = {'x': xcoords, 'y': ycoords}

    write_file(sys.argv[2], scales, coords, datapoints, gradient, intersect, fitcoords, gravitypoint)
