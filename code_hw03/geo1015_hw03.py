#-- geo1015_hw03.py
#-- Assignment 03 GEO1015/2018
#-- 2018-12-14

#------------------------------------------------------------------------------
# DO NOT MODIFY THIS FILE!!!
#------------------------------------------------------------------------------

import json 
import time
import rasterio

#-- all your code goes into my_code
import my_code_hw03


def main():
    
    #-- read the needed parameters from the file 'params.json' (must be in same folder)
    jparams = json.load(open('params2.json'))
    # jparams = json.load(open('params2.json'))

    start_time = time.time()

    #-- load in memory the input grid
    d = rasterio.open(jparams['input_file'])    

    #-- fetch the viewpoints
    viewpoints = []
    for i,each in enumerate(jparams['viewpoints']):
        vp = (jparams['viewpoints'][i]['xy'][0], jparams['viewpoints'][i]['xy'][1], jparams['viewpoints'][i]['height'])
        viewpoints.append(vp)

    #-- my code
    my_code_hw03.output_viewshed(d, viewpoints, jparams['maxdistance'], jparams['output_file'])

    print("--- %.3f seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()





