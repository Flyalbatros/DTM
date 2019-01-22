#-- my_code_hw04.py
#-- Assignment 04 GEO1015/2018
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 

import numpy as np
from scipy.linalg import eigh
from laspy import file
import scipy.spatial

import random
import math

def fit_plane(pts):
    """
    Fits a plane through a set of points using principal component analysis. 
     
    Input:
        pts:    the points to fit a plane through
        
    Output:
        (n,c)   a tuple with a point p that lies on the plane and the normalised normal vector n of the plane.
    """

    # shift points to mean
    mean = np.mean(pts, axis = 0)  
    pts -= mean
    # compute covariance matrix and eigenvalues and eignevectors
    cov = np.cov(pts, rowvar = False)
    evals, evecs = eigh(cov)
    #print("eigenvalues",evals)
    #print("eignevectors",evecs)
    # find smallest eigenvalue, the corresponging eigenvector is our normal n
    idx = np.argsort(evals)[::-1]
    #print("idx",idx)
    evecs = evecs[:,idx] #eigenvector
    evals = evals[idx]  #eigenvalue
    #print("evals",evals)
    #print("evecs",evecs)
    n = evecs[-1] #normal
    nor = evecs[:,-1]
    c = mean
    #print("n",n)
    #print("nor",nor)
    return c, nor

def ReadAndThinning(fn,thinning_factor):
    #returen a list of points after thinging 
    #pts [[x_coord,y_coord,z_coord].....]
    # [float,float,float]
    f = file.File(fn, mode='r')
    h = f.header
    #print("header\n")
    scale = h.scale
    offset = h.offset
    #print(h.scale,h.offset)
    #[0.01, 0.01, 0.01] [-0.0, -0.0, -0.0]

    noGroundPts = f.points[f.classification != 2]
    n = len(noGroundPts)
    pts =[]
    for i in list(range(0,n,thinning_factor)):
        p = noGroundPts[i]
        pts.append([p[0][0]*scale[0]+offset[0],p[0][1]*scale[1]+offset[1],p[0][2]*scale[2]+offset[2]])
    return pts

def write2PLY(fn,planes):
    #planes ndarray [[x,y,z,segment_id],......[]]
    print("start writePLY" + fn)
    f = open(fn, 'w')
    f.write("ply\n")
    f.write("format ascii 1.0\n")
    f.write("comment GEO1015 hw04\n")
    f.write("element vertex {}\n".format(len(planes)))
    f.write("property double x\n") 
    f.write("property double y\n")
    f.write("property double z\n")
    f.write("property int segment_id\n")
    f.write("end_header\n")
    np.savetxt(f, planes,fmt='%.2f %.2f %.2f %i', delimiter=' ')
    f.close()
    print("writePLY done")

def detect_planes(jparams):
    """
    !!! TO BE COMPLETED !!!
     
    Function that reads a LAS file, performs a region growing plane segmentation, and outputs the result to a PLY file.
     
    Input:
        jparams:  a dictionary with the paramameters:
            input-file:
                        the input LAS file
            output-file:
                        the output PLY file
            thinning-factor:        
                        thinning factor used to thin the input points prior to giving them to the segmentation 
                        algorithm. A factor of 1x means no thinning, 2x means remove 1/2 of the points, 
                        10x means remove 9/10 of the points, etc.
            minimal-segment-count:  
                        the minimal number of points in a segment.
            epsilon:                
                        distance threshold to be used during the region growing. It is the distance between 
                        a candidate point and the plane of the growing region.
            neighbourhood-type:     
                        specifies the type of neighbourhood search to be used, valid values are `knn` and `radius`.
            k:                      
                        size of the knn neighbourhood searches.
            r:                      
                        radius for he fixed-radius neighbourhood searches.
        
    Output:
        none (but output PLY file written to 'output-file')
    """  
    pts = ReadAndThinning(jparams['input-file'],jparams['thinning-factor'])
    '''
    unclassified_pts=[]
    for pt in pts:
        unclassified_pts.append([pt[0],pt[1],pt[2],0])
    ori_pts = np.array(unclassified_pts)
    write2PLY("thinning_notGround.ply",ori_pts)
    '''
    #n = len(pts)
    #print(n,int(n/50))
    pts = pts[:30000]
    #kd = scipy.spatial.KDTree([i[0:3] for i in pts])
    kd = scipy.spatial.KDTree(pts)

    #region growing
    #SeedInd = random.sample(list(range(0,len(pts))),  10*jparams['minimal-segment-count'])
    #SeedInd = list(range(0,len(pts),10*jparams['minimal-segment-count']))
    #print(SeedInd)
    remaining_pts = list(range(0,len(pts)))
    noPlane_pts = []
    out_plane = []
    plane_id = 0

    while(bool(remaining_pts)):
        current_plane = []
        current_plane_id =[]
        seed_stack = []

        randomSeed = int(random.sample(remaining_pts,1)[0])
        #randomSeed = 16349
        print("randomSeed -- ",randomSeed)
        seed_stack.append(randomSeed)
        current_plane.append(pts[randomSeed])
        current_plane_id.append(randomSeed)
        remaining_pts.remove(randomSeed)

        while(bool(seed_stack)):
            print("1",seed_stack)
            seed_id = seed_stack.pop()
            print("2",seed_stack)

            if jparams['neighbourhood-type'] == 'knn':
                distance, neighbourhoods_ind =  kd.query(pts[seed_id], k=jparams['k'])
            else:
                neighbourhoods_ind = kd.query_ball_point(pts[seed_id],jparams['r'])
            #print(neighbourhoods_ind)
            for p in neighbourhoods_ind:
                #print("1p",p)
                if p in remaining_pts :
                    #print("p in remaining_pts",p)
                    if len(current_plane)<3:
                        seed_stack.append(p)
                        current_plane.append(pts[p])
                        current_plane_id.append(p)
                        #print(p,type(p))
                        #print(remaining_pts[p-3:p+3])
                        remaining_pts.remove(p)

                    else:
                        print("start growing region")
                        testpt = pts[p]
                        current_plane_cen,current_plane_nor = fit_plane(current_plane)
                        v2 = np.array([testpt[0]-current_plane_cen[0],testpt[1]-current_plane_cen[1],testpt[2]-current_plane_cen[2]])
                        #Lv1 = np.sqrt(current_plane_nor.dot(current_plane_nor))
                        #because Lv1 == 1
                        #pt2plane = current_plane_nor.dot(v2)/Lv1
                        pt2plane = current_plane_nor.dot(v2)
                        if pt2plane  < jparams['epsilon']:
                            seed_stack.append(p)
                            current_plane.append(pts[p])
                            current_plane_id.append(p)
                            remaining_pts.remove(p)
            print("after growing -- seed_stack",len(seed_stack),seed_stack)

        if len(current_plane) > jparams['minimal-segment-count']:
            plane_id += 1
            print(plane_id,len(current_plane_id),current_plane_id)
            for pt in current_plane:
                out_plane.append([pt[0],pt[1],pt[2],plane_id])
        else:
            remaining_pts.extend(current_plane_id[1:])
            noPlane_pts.append(randomSeed)

    for ind in noPlane_pts:
        out_plane.append([pts[ind][0],pts[ind][1],pts[ind][2],-1])

    print("plane_id",plane_id)
 
    #planes = [[1.11,2.22,3.33,4],[1.23,2.34,3.23,4]]
    #planes = np.array([[1.11,2.22,3.33,4],[1.23,2.34,3.23,4]])
    planes = np.array(out_plane)
    #print(planes)
    write2PLY(jparams['output-file'],planes)
    

