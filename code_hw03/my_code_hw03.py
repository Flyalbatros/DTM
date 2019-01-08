#-- my_code_hw03.py
#-- Assignment 03 GEO1015/2018
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 

import sys
import math
import numpy
import rasterio
from rasterio import features



def output_viewshed(d, viewpoints, maxdistance, output_file):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster
     
    Input:
        d:            the input datasets (rasterio format)  
        viewpoints:   a list of the viewpoints (x, y, height)
        maxdistance:  max distance one can see
        output_file:  path of the file to write as output
        
    Output:
        none (but output GeoTIFF file written to 'output-file')
    """  
    #d --- <open DatasetReader name='data/tasmania/dem_01.tif' mode='r'>
    #viewpoints --- [(515709.0, 5263774.0, 2), (515900.0, 5263500.0, 2)]
    #maxdistance --- 2000
    
    # [this code can and should be removed/modified/reutilised]
    # [it's just there to help you]

    #-- numpy of input
    npi  = d.read(1)
    #print("npi",npi)
    '''
    [[750. 751. 750. ...  61.  62.  62.]
    [749. 751. 753. ...  60.  61.  61.]
    [748. 751. 754. ...  60.  60.  61.]
    ...
    [230. 236. 240. ...  11.  12.  12.]
    [234. 241. 247. ...  11.  13.  13.]
    [242. 248. 255. ...  12.  13.  14.]]
    '''

    PixelSizeX = d.transform[0]
    PixelSizeY = -d.transform[4]
    (nrow,ncol) = d.shape
    #print(PixelSizeX,PixelSizeY,shape)
    #33.365177999999986 33.36517799999934 (505, 550)
    
    maxPixels_x =  math.ceil(maxdistance/PixelSizeX)
    maxPixels_y = math.ceil(maxdistance/PixelSizeY)



    #-- fetch the 1st viewpoint
    v = viewpoints[0]
    #print(v)  (515709.0, 5263774.0, 2)
    #-- index of this point in the numpy raster
    vrow, vcol = d.index(v[0], v[1])
    #print(vrow,vcol)   340 320

    x_ind_min = vrow - maxPixels_x
    x_ind_max = vrow + maxPixels_x
    y_ind_min = vcol - maxPixels_y
    y_ind_max = vcol + maxPixels_y
    
    v1 = (x_ind_max, y_ind_max)
    v2 = (x_ind_max, y_ind_min)
    v3 = (x_ind_min, y_ind_min)
    v4 = (x_ind_min, y_ind_max)
    #print(v1,v2,v3,v4)
    #(400, 380) (400, 260) (280, 260) (280, 380)

    #-- the results of the viewshed in npvs, all values=0
    npvs = numpy.ones(d.shape, dtype=numpy.int8)*(3)
    ######################################################## check ! there is no 3 inside the R
    #print(npvs)
    #-- put that pixel with value 2
    npvs[vrow , vcol] = 2
    
    #creat r bounding box
    y = list(range(y_ind_min,y_ind_max+1))*2
    x = [x_ind_max]*(maxPixels_y*2)+[x_ind_min]*(maxPixels_y*2)
    x += list(range(x_ind_min+1,x_ind_max+1))*2
    y += [y_ind_max]*(maxPixels_x*2)+[y_ind_min]*(maxPixels_x*2)
    rBox = list(zip(x,y))[1:]
    #print(rBox)
    #[(400, 261), (400, 262), (400, 263),......]


    getOrderedIndList(d, 0, 0, (5,4))
    
    for item in rBox[0:1]:
        #print(item)
        getOrderedIndList(d, vrow, vcol, item)
        

    #for i in range(x_ind_min,x_ind_max):
    #    for j in range(y_ind_min,y_ind_max):
    #        central = 






    #-- write this to disk
    with rasterio.open(output_file, 'w', 
                       driver='GTiff', 
                       height=npi.shape[0],
                       width=npi.shape[1], 
                       count=1, 
                       dtype=rasterio.uint8,
                       crs=d.crs, 
                       transform=d.transform) as dst:
        dst.write(npvs.astype(rasterio.uint8), 1)

    print("Viewshed file written to '%s'" % output_file)

def getOrderedIndList(d, vrow, vcol, endcell):
    # d = rasterio dataset as above
    #a = (10, 10)
    #b = (100, 50)

    #-- create in-memory a simple GeoJSON LineString
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(d.xy(vcol,vrow))
    v["coordinates"].append(d.xy(endcell[0], endcell[1]))
    shapes = [(v, 1)]
    re = features.rasterize(shapes, 
                            out_shape=d.shape, 
                            # all_touched=True,
                            transform=d.transform)
    # re is a numpy with d.shape where the line is rasterised (values != 0)
    '''
    [[0 0 0 ... 0 0 0]
     [0 0 0 ... 0 0 0]
     [0 0 0 ... 0 0 0]
     ...
     [0 0 0 ... 0 0 0]
     [0 0 0 ... 0 0 0]
     [0 0 0 ... 0 0 0]]
    '''
    ansList = numpy.argwhere(re == 1)
    print(ansList)
    ansList.tolist()
    print(type(ansList),type(ansList[0]),type(ansList[0][1]))
    #<class 'numpy.ndarray'> <class 'numpy.ndarray'> <class 'numpy.int64'>
    '''
    if  endcell[1]>vrow :
        return reversed(ansList)
    else:
        return ansList
    '''



