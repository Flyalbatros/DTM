#-- my_code_hw01.py
#-- hw01 GEO1015/2018
#-- Pablo Ruben
#-- 4273818
#-- Qu Wang
#-- [YOUR STUDENT NUMBER] 


import math
import scipy.spatial
import numpy


def nn_interpolation(list_pts_3d, j_nn):
    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    raster_pts = compute_raster(min, max, j_nn['cellsize'])
    kd = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    d,i = kd.query([i[0:2] for i in list_pts_3d], k=1)
    #print ("distances",i)
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_nn:        the parameters of the input for "nn"
    Output:
        returns the value of the area
 
    """  
    print("=== Nearest neighbour interpolation ===")

    # print("cellsize:", j_nn['cellsize'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # d, i = kd.query(p, k=1)

    print("File written to", j_nn['output-file'])



def idw_interpolation(list_pts_3d, j_idw):
    print("=== IDW interpolation ===")
    #index the data
    kd = scipy.spatial.KDTree([i[0:2] for i in list_pts_3d])
    #compute raster
    min, max = compute_bbox([i[0:2] for i in list_pts_3d])
    raster_pts, line_length = compute_raster(min, max, j_idw['cellsize'])
    #calculate the idw values for each raster point
    print("radius:", j_idw['radius'])
    i = kd.query_ball_point(raster_pts, j_idw['radius'])
    #print(len(raster_pts), len(i))
    power = j_idw['power']
    counter = 0
    out_matrix = []
    outline = []
    for idx in range(0,len(raster_pts)):
        interpolation = 0
        weight_sum = 0
        for point_idx in i[idx]:
            raster_pt_x = raster_pts[idx][0]
            raster_pt_y = raster_pts[idx][1]
            interp_pt_x = list_pts_3d[point_idx][0]
            interp_pt_y = list_pts_3d[point_idx][1]
            distance = (((raster_pt_x-interp_pt_x)**2+(raster_pt_y-interp_pt_y)**2)**0.5)
            if distance == 0:
                interpolation = list_pts_3d[point_idx][2]
                weight_sum = 1
                break
            weight = distance**-power
            interpolation += list_pts_3d[point_idx][2]*weight
            weight_sum += weight
        interpolation = interpolation/weight_sum
        outline.append(interpolation)
        counter +=1
        if counter == line_length:
            counter = 0
            out_matrix.append(outline)
            outline = []
    print(out_matrix)


    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with IDW
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_idw:       the parameters of the input for "idw"
    Output:
        returns the value of the area
 
    """

    # print("cellsize:", j_idw['cellsize'])
    # print("radius:", j_idw['radius'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # i = kd.query_ball_point(p, radius)
    
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

    #-- example to construct the DT
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
    return (min_x[0],min_y[0]),(max_x[0], max_y[0])

def compute_raster(min, max, cellsize):
    out_points = []
    y_axis = range(int(min[1]), int(max[1])+int(cellsize), +int(cellsize))
    x_axis = range(int(min[0]), int(max[0])+int(cellsize), int(cellsize))
    for line in reversed(y_axis):
        for cell in x_axis:
            #print(cell,line)
            out_points.append((cell, line))
    return out_points, len(x_axis)