# -- my_code_hw01.py
# -- hw01 GEO1015/2018
# -- Pablo Ruben
# -- 4273818
# -- Qu Wang
# -- 4700686

import math
import scipy.spatial
import numpy as np
from numpy.linalg import inv


def nn_interpolation(list_pts_3d, j_nn):
    print("=== Nearest neighbour interpolation ===")
    cellsize = j_nn['cellsize']

    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    xll = min[0]
    yll = min[1]
    (ncols, nrows, rcpt) = compute_rcpt(min, max, cellsize)

    kd_tree = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    dd, ii = kd_tree.query(rcpt, k=1)

    rasterdata = np.ones((nrows, ncols))
    rasterdata = rasterdata * (-9999.0)

    ind_rcpt = 0
    for ind in ii:
        z = list_pts_3d[ind][2]
        i_row = int(ind_rcpt/ncols)
        j_col = ind_rcpt%ncols
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
    #print("radius:", j_idw['radius'])
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
            #print(interpolation, weight_sum)
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
    print("=== TIN interpolation ===")
    #clean the data
    list_pts_3d=clean_points(list_pts_3d)

    cellsize = j_tin['cellsize']
    #generate raster
    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    xll = min[0]
    yll = min[1]
    (ncols, nrows, rcpt) = compute_rcpt(min, max, cellsize)
    dt = scipy.spatial.Delaunay(list_pts_3d)
    #prepare output
    rasterdata = np.ones((nrows, ncols))
    rasterdata = rasterdata * (-9999.0)

    for cpt in rcpt:
        find_simplex()


    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html#scipy.spatial.Delaunay

    dt = scipy.spatial.Delaunay([])

    print("File written to", j_tin['output-file'])


def kriging_interpolation(list_pts_3d, j_kriging):
    print("=== Ordinary kriging interpolation ===")
    #clean the data
    list_pts_3d=clean_points(list_pts_3d)
    cellsize = j_kriging['cellsize']
    #generate raster
    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    xll = min[0]
    yll = min[1]
    (ncols, nrows, rcpt) = compute_rcpt(min, max, cellsize)
    #prepare output
    rasterdata = np.ones((nrows, ncols))
    rasterdata = rasterdata * (-9999.0)
    #kd-tree indexing
    kd_tree = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    #for each centerpoint
    r = j_kriging['radius']
    cpt_ind = 0
    for cpt in rcpt:
        #neighbours -- a list of the indices of the neighbors of x
        neighbours = kd_tree.query_ball_point(cpt, r) 
        #print("cpt",cpt)
        #print("len(neighbours)",len(neighbours),"neighbours",neighbours)
        if len(neighbours)>0:
            #build the a and d matrix
            A = np.ones((len(neighbours)+1, len(neighbours)+1))
            A[len(neighbours)][len(neighbours)] = 0
            d = np.ones((len(neighbours)+1, 1))
            for pt_ind1 in range(0,len(neighbours)):
                #traverse the list of neighbors
                x0 = cpt[0]
                y0 = cpt[1]
                x1 = list_pts_3d[neighbours[pt_ind1]][0]
                y1 = list_pts_3d[neighbours[pt_ind1]][1]
                distance = math.sqrt((x1-x0)**2+(y1-y0)**2)
                #print(distance)
                gamma_d = gaussian(distance)
                #write the distance to neighbor to d matrix
                d[pt_ind1] = gamma_d

                for pt_ind2 in range(pt_ind1,len(neighbours)):
                    #for all pairs of neighboring points
                    # calculate distance and use function to store in A matrix
                    x2 = list_pts_3d[neighbours[pt_ind2]][0]
                    y2 = list_pts_3d[neighbours[pt_ind2]][1]
                    distance = math.sqrt((x1-x2)**2+(y1-y2)**2)
                    gamma = gaussian(distance)
                    #print(gamma)
                    A[pt_ind1][pt_ind2]=gamma
                    A[pt_ind2][pt_ind1]=gamma
            #calculate weight by inverting matrix and scalar product
            w = np.dot(inv(A),d)
            #calculate the value of the point using the weights
            z = 0
            for i in range(0,len(neighbours)):
                z+=w[i]*list_pts_3d[neighbours[i]][2]
            #print("z",z)
            #store the result
            i_row = int(cpt_ind/ncols)
            j_col = cpt_ind%ncols
            rasterdata[i_row][j_col] = z
        cpt_ind+=1

    fname = j_kriging['output-file']
    writeASC(fname, ncols, nrows, xll, yll, cellsize, rasterdata)
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
    for y_row in reversed(range(0, nrows)):
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

def gaussian(distance, sill=1300, range=270, nugget=0):
#def gaussian(distance, sill=1450, range=300, nugget=0):
    #return sill*(1-np.exp(((-3*distance)**2)/(range^2)))+nugget
    return sill*(1-math.exp(-(3*distance)**2/range**2))+nugget


def clean_points(points_list):
    clean_points_list = []
    for point1 in points_list:
        repeated = False
        for point2 in clean_points_list:
            if point1[0] == point2[0] and point1[1] == point2[1]:
                repeated = True
        if repeated == False:
            clean_points_list.append(point1)
        else:
            print("Repeated point: " + str(point1[0]) + " " + str(point1[1]))
    return np.array(clean_points_list)
