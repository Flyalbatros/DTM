#-- geo1015_hw02.py
#-- hw02 GEO1015/2018
#-- 2018-12-03

#------------------------------------------------------------------------------
# DO NOT MODIFY THIS FILE!!!
#------------------------------------------------------------------------------



import json, sys

from scipy.spatial import Delaunay
import numpy as np

from my_code_hw02 import read_pts_from_grid, simplify_by_refinement, compute_differences

def write_obj(pts, obj_filename):
    print("=== Writing {} ===".format(obj_filename))
    f_out = open(obj_filename, 'w')
    for p in pts:
        f_out.write("v {:.2f} {:.2f} {:.2f}\n".format(p[0], p[1], p[2]))
    dt = Delaunay(pts[:,:2])
    for simplex in dt.simplices:
        f_out.write("f {} {} {}\n".format(simplex[0]+1, simplex[1]+1, simplex[2]+1))
    f_out.close()

def main():
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()

    # step 1: convert the grid to a set of points
    pts = read_pts_from_grid(jparams)

    # step 2: select the important points from pts using the TIN refinement algorithm
    pts_important = simplify_by_refinement(pts, jparams)
    # write the simplified TIN to OBJ:
    if pts_important is not None:
        write_obj(pts_important, jparams['output-file-tin'])
    else:
        print("No OBJ written since no important points.")

    # step 3: compare the simplified TIN with the original input grid and write the per-pixel differences to a new grid
    compute_differences(pts_important, jparams)

if __name__ == "__main__":
    main()