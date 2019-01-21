#-- my_code_hw04.py
#-- Assignment 04 GEO1015/2018
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 
#-- [YOUR NAME] 
#-- [YOUR STUDENT NUMBER] 

import numpy as np
import scipy.spatial
from scipy.linalg import eigh
from laspy import file
import random



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
    # find smallest eigenvalue, the corresponging eigenvector is our normal n
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    evals = evals[idx]

    n = evecs[-1]
    c = mean

    return n,c

def valid_candidate(plane_pt, normal_v, query_pt, epsilon):
    plane_query_v = np.array([plane_pt[0]-query_pt[0],plane_pt[1]-query_pt[1],plane_pt[2]-query_pt[2]])
    return np.dot(plane_query_v,normal_v)<epsilon

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
    #filter out the non ground points
    noGroundPts = f.points[f.classification != 2]
    n = len(noGroundPts)
    pts =[]
    #thin and reproject the points
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
    tree = scipy.spatial.cKDTree(pts)
    #print(pts[0:10])
    #print(len(pts))

    #region growing
    #SeedInd = random.sample(list(range(0,len(pts))),  10*jparams['minimal-segment-count'])
    SeedInd = list(range(0,len(pts),10*jparams['minimal-segment-count']))
    #print(SeedInd)

    #create a dictionarry to store the visited points and the region they belong to
    region_dict = {}
    for reg_index in range(0,len(SeedInd)):
        curr_pt = pts[SeedInd[reg_index]]
        curr_region = {SeedInd[reg_index]:reg_index}
        #query the neighbors of the seed
        neighb = neighb_querier(jparams, curr_pt, tree)
        #determine if seed is a good one
        seed_neighs = []
        for n_ind in neighb:
            n_pt = pts[n_ind]
            seed_neighs.append(n_pt)
        normal_v, plane_pt = fit_plane(seed_neighs)
        #print("first", neighb)
        #let's identify the points of the first expansion
        first_expansion = expand_region(neighb, plane_pt, normal_v, pts, jparams)
        if first_expansion != -1:
            filtered = neighb_filter(first_expansion, curr_region, [], region_dict)
            neighb = filtered
        expansion_ind = -99
        #let's perform iterative region growing until no points can be added anymore
        while expansion_ind !=-1:
            #print("second", neighb)
            expansion_ind = expand_region(neighb, plane_pt, normal_v, pts, jparams)
            #only if the region can be grown...
            if expansion_ind != -1:
                #add the new points to the current regions
                expansion_pts = []
                for ind in expansion_ind:
                    expansion_pts.append(pts[ind])
                    if (ind in curr_region) == False:
                        curr_region[ind] = reg_index
                #get the new neighbors for the next step, but avoid including existing points
                neighb = []
                for pt in expansion_pts:
                    neighb_inds = neighb_querier(jparams, pt, tree)
                    filtered = neighb_filter(neighb_inds, curr_region, neighb, region_dict)
                    neighb += filtered
                # if neighb == []:
                #     expansion_ind = -1
                #print("third", neighb, curr_region)
        if len(curr_region)>jparams["minimal-segment-count"]:
            region_dict = {**curr_region, **region_dict}
        print(region_dict)

def neighb_filter(neighb_inds, curr_region, neighb, region_dict):
    outlist = []
    for neighb_ind in neighb_inds:
        if (neighb_ind in curr_region) == False and (neighb_ind in neighb) == False and (neighb_ind in region_dict) == False:
            outlist.append(neighb_ind)
    return outlist

def expand_region(neighb, plane_pt, normal_v, pts, jparams):
    expansion = []
    #print(neighb)
    for query_pt_ind in neighb:
        #print("iterating")
        if valid_candidate(plane_pt, normal_v, pts[query_pt_ind], jparams['epsilon'])==True:
            expansion.append(query_pt_ind)
            #print("added point")
        #else:
            #print("point too far", query_pt_ind)
    if len(expansion) == 0:
        #print("no expansion possible, trying next seed/neighbor point")
        return -1
    else:
        return expansion

def neighb_querier(jparams, curr_pt, tree):
    #print("current", curr_pt)
    if type(curr_pt[0]) == type([]):
        print("error: more than one point provided for neighbour query")
        return "error"
    if jparams['neighbourhood-type'] == 'knn':
        dist, neighb = tree.query(curr_pt,jparams['k'])
    if jparams['neighbourhood-type'] == 'radius':
        neighb = tree.query_ball_point(curr_pt,jparams['r'])
    return neighb










    #planes = [[1.11,2.22,3.33,4],[1.23,2.34,3.23,4]]
    planes = np.array([[1.11,2.22,3.33,4],[1.23,2.34,3.23,4]])
    print(planes)
    write2PLY(jparams['output-file'],planes)


