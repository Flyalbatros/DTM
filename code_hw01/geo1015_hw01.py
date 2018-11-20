#-- geo1015_hw01.py
#-- hw01 GEO1015/2018
#-- Hugo Ledoux <h.ledoux@tudelft.nl>
#-- 2018-11-07

#------------------------------------------------------------------------------
# DO NOT MODIFY THIS FILE!!!
#------------------------------------------------------------------------------

import sys
import math
import csv
import random
import json 

#-- *all* your code goes into 'my_code_hw01'
import my_code_hw01


def main():
    #-- read the needed parameters from the file 'params.json' (must be in same folder)
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()
    #-- store the input 3D points in list
    list_pts_3d = []
    with open(jparams['input-file']) as csvfile:
        r = csv.reader(csvfile, delimiter=' ')
        header = next(r)
        for line in r:
            p = list(map(float, line)) #-- convert each str to a float
            assert(len(p) == 3)
            list_pts_3d.append(p)
    #-- interpolations if in the params
    if 'nn' in jparams:
        my_code_hw01.nn_interpolation(list_pts_3d, jparams['nn'])
    if 'idw' in jparams:
        my_code_hw01.idw_interpolation(list_pts_3d, jparams['idw'])
    if 'tin' in jparams:
        my_code_hw01.tin_interpolation(list_pts_3d, jparams['tin'])
    if 'kriging' in jparams:
        my_code_hw01.kriging_interpolation(list_pts_3d, jparams['kriging'])

if __name__ == '__main__':
    main()





