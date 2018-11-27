# -- my_code_hw01.py
# -- hw01 GEO1015/2018
# -- Pablo Ruben
# -- 4273818
# -- Qu Wang
# -- 4700686

import math
import scipy.spatial
import numpy as np


def nn_interpolation(list_pts_3d, j_nn):
    print("=== Nearest neighbour interpolation ===")
    cellsize = j_nn['cellsize']

    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    xll = min[0]
    yll = min[1]

    (ncols, nrows, rcpt) = compute_rcpt(min, max, cellsize)

    # xxyy = [i[0:2] for i in list_pts_3d]
    # print("xxyy",xxyy[0:10])
    kd_tree = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    # print("kd_tree",kd_tree.data[0:10])

    dd, ii = kd_tree.query(rcpt, k=1)

    rasterdata = np.ones((nrows, ncols))
    rasterdata = rasterdata * (-9999.0)

    # d = pt3d2dic(list_pts_3d)
    # print("nrows*ncols",nrows*ncols)
    # print("len(rcpt)",len(rcpt))
    # print("len(ii)",len(ii))

    ind_rcpt = 0
    for ind in ii:
        z = list_pts_3d[ind][2]
        # print ("z",z)
        i_row = math.floor((rcpt[ind_rcpt][1] - yll) / cellsize)
        j_col = math.floor((rcpt[ind_rcpt][0] - xll) / cellsize)
        rasterdata[i_row][j_col] = z
        ind_rcpt += 1

    fname = j_nn['output-file']
    writeASC(fname, ncols, nrows, xll, yll, cellsize, rasterdata)
    print("File written to", j_nn['output-file'])

def idw_interpolation(list_pts_3d, j_idw):
    print("=== IDW interpolation ===")
    # index the data
    kd = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    # compute raster
    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    (ncols, nrows, raster_pts) = compute_rcpt(min, max, j_idw['cellsize'])
    # calculate the idw values for each raster point
    print("radius:", j_idw['radius'])
    #print(raster_pts)
    i = kd.query_ball_point(raster_pts, j_idw['radius'])
    # print(len(raster_pts), len(i))
    power = j_idw['power']
    column_counter = 0
    line_counter = 0
    out_matrix = np.ones((nrows, ncols))
    out_matrix = out_matrix * (-9999.0)
    outline = []
    for idx in range(0, len(raster_pts)):
        interpolation = 0
        weight_sum = 0
        if i[idx] == []:
            break
        for point_idx in i[idx]:
            raster_pt_x = raster_pts[idx][0]
            raster_pt_y = raster_pts[idx][1]
            interp_pt_x = list_pts_3d[point_idx][0]
            interp_pt_y = list_pts_3d[point_idx][1]
            distance = (((raster_pt_x - interp_pt_x) ** 2 + (raster_pt_y - interp_pt_y) ** 2) ** 0.5)
            if distance == 0:
                interpolation = list_pts_3d[point_idx][2]
                weight_sum = 1
                break
            weight = distance ** -power
            interpolation += list_pts_3d[point_idx][2] * weight
            weight_sum += weight
            print(interpolation, weight_sum)
        interpolation = interpolation / weight_sum
        out_matrix[line_counter][column_counter]=interpolation
        column_counter += 1
        if column_counter == ncols:
            column_counter = 0
            outline = []
            line_counter +=1
    fname = j_idw['output-file']
    writeASC(fname, ncols, nrows, min[0], min[1], j_idw['cellsize'], out_matrix)
    print("File written to", j_idw['output-file'])

def tin_interpolation(list_pts_3d, j_tin):
    """
    !!! TO BE COMPLETED !!!

    Function that writes the output raster with linear in TIN interpolation

    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_tin:       the parameters of the input for "tin"
    Output:
        returns the value of the area

    """
    print("=== TIN interpolation ===")

    # -- example to construct the DT
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html#scipy.spatial.Delaunay
    # dt = scipy.spatial.Delaunay([])

    print("File written to", j_tin['output-file'])


def kriging_interpolation(list_pts_3d, j_kriging):
    """
    !!! TO BE COMPLETED !!!

    Function that writes the output raster with ordinary kriging interpolation

    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_kriging:       the parameters of the input for "kriging"
    Output:
        returns the value of the area

    """
    print("=== Ordinary kriging interpolation ===")

    print("File written to", j_kriging['output-file'])


###additional functions###

def compute_bbox(list_pts):
    min_x = min([i[0:1] for i in list_pts])
    min_y = min([i[1:2] for i in list_pts])
    max_x = max([i[0:1] for i in list_pts])
    max_y = max([i[1:2] for i in list_pts])
    # print((min_x[0],min_y[0]),(max_x[0], max_y[0]))
    return (min_x[0], min_y[0]), (max_x[0], max_y[0])

def compute_rcpt(min, max, cellsize):
    min_x = min[0]
    min_y = min[1]
    max_x = max[0]
    max_y = max[1]
    ncols = math.ceil((max_x - min_x) / cellsize)  # x
    nrows = math.ceil((max_y - min_y) / cellsize)  # y
    rcpt = []
    for y_row in range(0, nrows):
        y_coor = min_y + (y_row + 0.5) * cellsize
        for x_col in range(0, ncols):
            x_coor = min_x + (x_col + 0.5) * cellsize
            central_pt = [x_coor, y_coor]
            rcpt.append(central_pt)
    return (ncols, nrows, rcpt)

def writeASC(fname, ncols, nrows, xll, yll, cellsize, rasterdata):
    print("start writeASC" + fname)
    f = open(fname, 'w')
    f.write("NCOLS {}\n".format(ncols))
    f.write("NROWS {}\n".format(nrows))
    f.write("XLLCORNER {}\n".format(xll))
    f.write("YLLCORNER {}\n".format(yll))
    f.write("CELLSIZE {}\n".format(cellsize))
    f.write("NODATA_VALUE -9999\n")
    np.savetxt(f, rasterdata, fmt='%.1f', delimiter=' ')
    f.close()
    print("writeASC done")
