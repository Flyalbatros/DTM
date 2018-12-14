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
    transform_matrix = input_data.transform
    print(transform_matrix)
    lower_left_corner = transform_matrix * (0,input_data.height)
    for row in range(0,nrows):
        for col in range(0,ncols):
            print(transform_matrix * (row,col),
            z=raw_data[0][row][col])
            pass
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

    # Remember: the vertices of the initial TIN should not be returned



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
    