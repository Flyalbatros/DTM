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
    #-- numpy of input
    npi  = d.read(1)
    # get the pixel sizes
    PixelSizeX = d.transform[0]
    PixelSizeY = -d.transform[4]
    # translate the range into pixels
    maxPixels_x =  math.ceil(maxdistance/PixelSizeX)
    maxPixels_y = math.ceil(maxdistance/PixelSizeY)
    #write 3 (outside view range, default value) on the entire output dataset
    npvs = numpy.ones(d.shape, dtype=numpy.int8) * (3)
    #for each viewpoint...
    for v in viewpoints:
        #determine the position in the raster and register it (value 2 in output raster)
        vrow, vcol = d.index(v[0], v[1])
        npvs[vrow, vcol] = 2
        #and the absolute altitude of the viewpoint...
        alt_vp = npi[(vrow, vcol)] + v[2]
        #determine the boundaries of the box circumscribed to which the circle is circumscribed
        x_ind_min = vrow - maxPixels_x
        x_ind_max = vrow + maxPixels_x
        y_ind_min = vcol - maxPixels_y
        y_ind_max = vcol + maxPixels_y
        #creat a squared bounding box using previous boundary values
        y = list(range(y_ind_min,y_ind_max+1))*2 #make sure to include the max in the range...
        x = [x_ind_max]*(maxPixels_y*2+1)+[x_ind_min]*(maxPixels_y*2+1)
        x += list(range(x_ind_min,x_ind_max+1))*2
        y += [y_ind_max]*(maxPixels_x*2+1)+[y_ind_min]*(maxPixels_x*2+1)
        rBox = list(zip(x,y))[1:] # this is the box
        #just save the square of the max distance so that it doesn't need to be computed on the fly
        sq_maxdistance = maxdistance**2
        #let's transform the box into a circle by shrinkening it
        cBox = []
        for index in range(0,len(rBox)):
              dist_pt = rBox[index]
              #calculate distance between viewpoint and border point on the box
              sq_dist = ((dist_pt[0]-vrow)*PixelSizeY)**2+((dist_pt[1]-vcol)*PixelSizeX)**2
              #if the distance is too long, shrinken the ray
              if sq_dist>sq_maxdistance:
                  ratio=(sq_maxdistance)**0.5/(sq_dist)**0.5
                  new_coords = (math.ceil(vrow+(dist_pt[0]-vrow)*ratio), math.ceil(vcol+(dist_pt[1]-vcol)*ratio))
              #if distance is not too high, the coordinates stay the same
              else:
                  new_coords = (dist_pt[0],dist_pt[1])
              #now let's write the circle boundary points and filter out any duplicates
              if not new_coords in cBox:
                  cBox.append(new_coords)
        #finally use Bresenheim to get the paths of the rays leading to the circle bounding points
        for bound_pt in cBox:
            path = Bresenham_with_rasterio(d, (vrow, vcol), bound_pt)
            #set the default values for the biggest altitude and tangent observed in the ray
            max_alt = -99
            max_tan = -99
            for point in path[1:]:
                #get the altitude of the point
                alt = npi[point[0],point[1]]
                #if the altitude is high enough...
                if alt>max_alt:
                    dist = calCosAng(path[0],path[-1],point)#get the distance in number of cells of the projected center point on the line
                    if alt>=alt_vp: # only change the max needed altitude if the obstructing pixel is higher than the original viewpoint
                        max_alt = alt
                    rel_alt = alt - alt_vp #calculate the relative altitude difference with regard to viewpoint altitude
                    if rel_alt>0: #if the viewpoint is lower
                        tan = numpy.arctan(rel_alt/(dist))
                    elif rel_alt<0: #if the viewpoint is higher
                        tan = numpy.arctan(-rel_alt/(dist))*(-1)
                    elif rel_alt==0: #if the viewpoint is at same height
                        tan=0
                    if tan>max_tan and npvs[point[0], point[1]] != 2: #if the current tangent is bigger than the biggest so far = visible point & make sure it is no viewpoint
                        max_tan=tan #update the value
                        npvs[point[0], point[1]] = 1 # mark point as visible
                    elif npvs[point[0], point[1]] == 3: #only if point is not visible yet
                        npvs[point[0], point[1]] = 0 #mark point as unvisible
                elif npvs[point[0], point[1]] == 3:
                    #if the altitude of the next point is lower or equal to the altitude of the highest previous point location above viewpoint height, no need to check
                    npvs[point[0], point[1]] = 0

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

def Bresenham_with_rasterio(raster, start, end):
    d = raster
    a = start #format is (row,col)
    b = end #format is (row, col)
    #-- create in-memory a simple GeoJSON LineString
    v = {}
    v["type"] = "LineString"
    v["coordinates"] = []
    v["coordinates"].append(d.xy(a[0], a[1]))
    v["coordinates"].append(d.xy(b[0], b[1]))
    shapes = [(v, 1)]
    re = features.rasterize(shapes,
                            out_shape=d.shape,
                            all_touched=True,
                            transform=d.transform)
    out = numpy.argwhere(re==1)
    outlist = []
    for el in out:
        outlist.append(tuple(el))
    #depending on the orientation of the line, sort cells in right order
    if a[0]>b[0] and a[1]<=b[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], x[1]))
        #print('a')
    elif a[0]>b[0] and a[1]>=b[1]:
        outlist = sorted(outlist, key=lambda x: (-x[0], -x[1]))
        #print('b')
    elif a[0]<=b[0] and a[1]>b[1]:
        outlist = sorted(outlist, key=lambda x: (x[0], -x[1]))
    return outlist


def calCosAng(start,end,current_pt):
    #extract the two vectors
    v1 = numpy.array([end[1]-start[1],end[0]-start[0]])
    v2 = numpy.array([current_pt[1]-start[1],current_pt[0]-start[0]])
    #get the length of the vector describing the line
    Lv1 = numpy.sqrt(v1.dot(v1))
    #and finally get the projected distance of current_pt on the line
    dist = v1.dot(v2)/Lv1
    return dist
