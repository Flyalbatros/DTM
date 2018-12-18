#-- my_code_hw02.py
#-- hw02 GEO1015/2018
#-- [YOUR NAME]
#-- [YOUR STUDENT NUMBER] 
#-- [YOUR NAME]
#-- [YOUR STUDENT NUMBER] 


import scipy.spatial
import numpy as np
import rasterio


def read_pts_from_grid(jparams):
    """
    !!! TO BE COMPLETED !!!
     
    Function that reads a grid in .tif format and retrieves the pixels as a list of (x,y,z) points shifted to the origin
     
    Input from jparams:
        input-file:  a string containing the path to a grid file with elevations
    Returns:
        a numpy array where each row is one (x,y,z) point. Each pixel from the grid gives one point (except for no-data pixels, these should be skipped).
    """
    print("=== Reading points from grid ===")
    input_data = rasterio.open(jparams["input-file"])
    nodata_value = input_data.nodata
    raw_data = input_data.read()
    ncols = input_data.width
    nrows = input_data.height
    shape = input_data.shape
    PixelSizeX = input_data.transform[0]
    PixelSizeY = -input_data.transform[4]
    #transform_matrix = input_data.transform
    #print(raw_data[0])
    #lower_left_corner = transform_matrix * (0,input_data.height)
    # outlist = []
    # row_cnt = 0
    # for rev_row in reversed(range(0,nrows)):
    #     for col in range(0,ncols):
    #         z=raw_data[0][row_cnt][col]
    #         if z!=nodata_value:
    #             outlist.append([(rev_row+0.5)*PixelSizeY, (col+0.5)*PixelSizeX, z])
    #     row_cnt+=1
    # print(len(raw_data[0][1]))
    # print(len(np.array(generate_raster_points(nrows, ncols, raw_data, nodata_value, PixelSizeX, PixelSizeY))))
    return np.array(generate_raster_points(nrows, ncols, raw_data, nodata_value, PixelSizeX, PixelSizeY))
    # Tip: the most efficient implementation of this function does not use any loops. Use numpy functions instead.



def simplify_by_refinement(pts, jparams):
    """
    !!! TO BE COMPLETED !!!
     
    Function that takes a list of points and constructs a TIN with as few points as possible, while still satisfying the error-threshold. 

    This should be an implemented as a TIN refinement algorithm using greedy insertion. As importance measure the vertical error should be used.
     
    Input:
        pts:                a numpy array with on each row one (x,y,z) point
        from jparams:
            error-threshold:    a float specifying the maximum allowable vertical error in the TIN
    Returns:
        a numpy array that is a subset of pts and contains the most important points with respect to the error-threshold
    """
    print("=== TIN simplification ===")
    bbox_size = 3 #variable for bounding box size
    y_max = max(pts[:,1])
    x_max = max(pts[:,0])
    y_min = min(pts[:,1])
    x_min = min(pts[:,0])
    y_delta = y_max-y_min
    x_delta = x_max-x_min
    y_max += y_delta*0.5*(bbox_size-1)
    y_min -= y_delta*0.5*(bbox_size-1)
    x_max += x_delta*0.5*(bbox_size-1)
    x_min -= x_delta*0.5*(bbox_size-1)
    z_avg = sum(pts[:,2])/len(pts[:,2])
    dt_vertices = np.array([[x_min,y_min,z_avg], [x_max, y_min,z_avg], [x_max, y_max,z_avg], [x_min, y_max,z_avg]])
    #print(dt_vertices)
    dt_2d = scipy.spatial.Delaunay([i[0:2] for i in dt_vertices])
    error_track = 0
    highest_diff = np.inf
    while highest_diff>jparams["error-threshold"] and error_track==0:
        diff_list = []
        for point in pts:
            triangle_idx = dt_2d.find_simplex(point[0:2])
            #print(triangle_idx)
            if triangle_idx == -1:
                print("!!! error creating the bounding box !!!")
                error_track = 1
                break
            else: #calculate the difference between the existing TIN and the actual z value of the point
                interpolation = TIN_interpolator(dt_vertices, dt_2d, triangle_idx,  point)
                diff_list.append(abs(point[2]-interpolation))
        #update values and triangulation
        highest_diff = max(diff_list)
        if highest_diff>jparams["error-threshold"]:
            max_idx = diff_list.index(max(diff_list))
            dt_vertices = np.append(dt_vertices,[pts[max_idx]], axis=0)
            dt_2d = scipy.spatial.Delaunay([i[0:2] for i in dt_vertices])
    #print("%.32f" %highest_diff)
    #print(max(diff_list), min(diff_list))
    if len(dt_vertices)>4:
        return dt_vertices[4:len(dt_vertices)] # Remember: the vertices of the initial TIN should not be returned
    else:
        return None


def compute_differences(pts_important, jparams):
    """
    !!! TO BE COMPLETED !!!
     
    Function that computes the elevation differences between the input grid and the Delaunay triangulation that is constructed from pts_important. The differences are computed for each pixel of the input grid by subtracting the grid elevation from the TIN elevation. The output is a new grid that stores these differences as float32 and has the same width, height, transform, crs and nodata value as the input grid.

    Input:
        pts_important:          numpy array with the vertices of the simplified TIN
        from jparams:
            input-file:                 string that specifies the input grid
            output-file-differences:    string that specifies where to write the output grid file with the differences
    """
    print("=== Computing differences ===")
    input_data = rasterio.open(jparams["input-file"])
    out_profile = input_data.profile
    out_profile['dtype'] = 'float32'
    raw_data = input_data.read()
    PixelSizeX = input_data.transform[0]
    PixelSizeY = -input_data.transform[4]
    ###
    nodata_value = input_data.nodata
    ncols = input_data.width
    nrows = input_data.height
    shape = input_data.shape
    ###
    raster_pts = np.array(generate_raster_points(nrows, ncols, raw_data, nodata_value, PixelSizeX, PixelSizeY,0))
    ### generate the simplified TIN
    dt_2d = scipy.spatial.Delaunay([i[0:2] for i in pts_important])
    ###now let's compare them
    outlist = []
    linelist = []
    # print(ncols,nrows)
    # print(shape)
    # print(len(raster_pts))
    # print(len(raw_data[0][1]))
    col_counter = 0
    row_counter = 0
    for point in raster_pts:
        if point[2] == nodata_value:
            linelist.append(nodata_value)
        else:
            triangle_idx = dt_2d.find_simplex(point[0:2])
            if triangle_idx == -1:
                print("!!! WARNING: point outside convex hull of simplified dataset !!!")
                linelist.append(nodata_value)
            else:
                interpolation = TIN_interpolator(pts_important, dt_2d, triangle_idx, point)
                linelist.append(point[2] - interpolation)
        #index counters
        col_counter +=1
        if col_counter == ncols:
            col_counter = 0
            outlist.append(linelist)
            linelist = []
    #print(diff_raster)
    #let's write the output file reusing the settings of the input file
    outputter = rasterio.open(jparams["output-file-differences"], 'w', **out_profile)
    outputter.write(np.array([outlist]).astype(rasterio.float32))

def TIN_interpolator(dt_vertices, dt_2d, triangle_idx,  point):
    tri_vertices = dt_2d.simplices[triangle_idx]
    vertex_1 = dt_vertices[tri_vertices[0]]
    vertex_2 = dt_vertices[tri_vertices[1]]
    vertex_3 = dt_vertices[tri_vertices[2]]
    p = np.array([point[0:2]])  # put the point into a numpy array
    b = dt_2d.transform[triangle_idx, :2].dot(np.transpose(p - dt_2d.transform[triangle_idx, 2]))  # calculate barycentric coordinates 1/2
    weights = np.c_[np.transpose(b), 1 - b.sum(axis=0)][0]  # calculate barycentric coordinates 2/2
    return vertex_1[2] * weights[0] + vertex_2[2] * weights[1] + vertex_3[2] * weights[2]

def generate_raster_points(nrows, ncols, raw_data, nodata_value, PixelSizeX, PixelSizeY, nodata_filter=1):
    outlist = []
    row_cnt = 0
    for rev_row in reversed(range(0, nrows)):
        for col in range(0, ncols):
            z = raw_data[0][row_cnt][col]
            if z != nodata_value and nodata_filter==1:
                outlist.append([(col + 0.5) * PixelSizeX, (rev_row + 0.5) * PixelSizeY, z])
            if nodata_filter==0:
                outlist.append([(col + 0.5) * PixelSizeX, (rev_row + 0.5) * PixelSizeY, z])
        row_cnt += 1
    #print(row_cnt)
    return outlist