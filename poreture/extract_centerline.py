# -*- coding: utf-8 -*-

"""
This module implements centerline extraction and skeletonization methods.

The methods work on 3D binary volumes composed of
    0: background
    1: object to be skeletonized/centerline extracted

Literature references for implemented methods
    kline_vessel - [Kline et al. ABME 2010]
    kline_pore - [Kline et al. J Porous Mat 2011]
"""

# Imports
import os
import sys
import time
import nibabel as nib
import numpy as np
import skfmm # pip install scikit-fmm
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

def find_3D_object_voxel_list(vol):
    """This function creates a centerline from the segmented volume (vol)

    Inputs:
        Required:
            vol: 3D binary volume where -
                0: background
                1: object to be skeletonized/centerline extracted
    Returns:
        nz_x, nz_y, nz_z
            lists of non-zero x, y, z co-ordinates

    Dependencies:
        numpy
    """
    nz_coords = np.nonzero(vol)
    nz_x, nz_y, nz_z = [nz_coords[i].tolist() for i in range(3)]
    return nz_x, nz_y, nz_z


def find_terminal_end_points(vol):

    Place_holder = 1

def detect_local_maxima(vol):
    """
    Takes a 3D volume and detects the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define a 26-connected neighborhood
    neighborhood = morphology.generate_binary_structure(3,3) # first is dimension, next is relative connectivity

    # apply the local maximum filter; all locations of maximum value 
    # in their neighborhood are set to 1
    local_max = (filters.maximum_filter(vol, footprint=neighborhood)==vol)

    # Remove background
    local_max[vol==0] = 0

    # Find endpoint indici
    [xOrig,yOrig,zOrig] = np.shape(vol)
    x = []
    y = []
    z = []
    for i in range(0,xOrig):
        for j in range(0,yOrig):
            for k in range(0,zOrig):
                if local_max[i,j,k] > 0:
                    x.append(i)
                    y.append(j)
                    z.append(k)

    return x, y, z

def kline_vessel(vol, startID, **kwargs):
    """This function creates a centerline from the segmented volume (vol)

    Inputs:
        Required:
            vol: 3D binary volume where -
                0: background
                1: object to be skeletonized/centerline extracted
            startID: index of root, as in [x,y,z] location
    
        Optional (kwargs):
            dist_map_weight
            cluster_graph_weight
            min_branch_length
            min_branch_to_root

    Returns:
        extracted_centerline

    dependencies:
        numpy
        scikit-fmm (import name: skfmm)
        scipy

    """
    
    # Make sure object is equal to 1, without specifying dtype, would be logical
    B2 = np.array(vol.copy() > 0, dtype = 'int8')

    # Set defaults
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            print "%s == %s" %(key,value)
            if key=='dist_map_weight':
                dmw = value
            if key=='cluster_graph_weight':
                cgw = value
            if key=='min_branch_length':
                mbl = value
            if key=='min_branch_to_root':
                mbtr = value
    if 'dmw' not in locals():
        dmw = 6
    if 'cgw' not in locals():
        cgw = np.sum(vol)/20
    if 'mbl' not in locals():
        mbl = 5
    if 'mbtr' not in locals():
        mbtr = 10
    print "dmw = %s" %(dmw)
    print "cgw = %s" %(cgw)

    # Remember original volume size
    [xOrig,yOrig,zOrig] = np.shape(B2)

    # Find 3D coordinates of volume

    x3, y3, z3 = find_3D_object_voxel_list(B2)

    # Limit volume size
    B2 = B2[np.min(x3):np.max(x3)+1,np.min(y3):np.max(y3)+1,np.min(z3):np.max(z3)+1]

    # Setup starting index list and correct for change in volume size
    sx = startID[0] - np.min(x3) 
    sy = startID[1] - np.min(y3) 
    sz = startID[2] - np.min(z3) 

    # New volume size (bounding box)
    [x_si,y_si,z_si] = np.shape(B2)

    sys.stdout.flush() 
    time.sleep(1)       
    # Perform first fast march to determine endpoints
    # works on binary speed function
    phi = B2.copy()
    constant_speed = B2.copy()
    phi[sx,sy,sz] = -1
    #constant_speed = np.ones(shape = (np.shape(phi)))
    mask = B2<1
    phi = np.ma.MaskedArray(phi, mask)
    binary_travel_time = skfmm.travel_time(phi, constant_speed)

    # Fill in masked values and set to zero
    binary_travel_time = binary_travel_time.filled()
    binary_travel_time[binary_travel_time==1.e20] = 0
    print "minimum of binary travel time is %s" % np.min(binary_travel_time)

    sys.stdout.flush() 
    time.sleep(1)       
    # Normalize and apply cluster graph weighting (cluster graph weighting doesn't seem to be doing much, perhaps a better FMM implementation???)
    # Find endpoints
    hold_binary_travel_time = binary_travel_time.copy()
    print "number of non-zero elements is %s" % (np.sum(B2))
    [endx, endy, endz] = detect_local_maxima(hold_binary_travel_time)
    print "number of local maxima was %s" % (len(endx))

    sys.stdout.flush() 
    time.sleep(1)       
    # Now perform second FMM, to create field for gradient descent
    dMap = morphology.distance_transform_edt(constant_speed) #distance map finds distance from 1's, to nearest 0.
    weighted_speed = dMap ** dmw
    weighted_travel_time = skfmm.travel_time(phi, weighted_speed)
    weighted_travel_time = weighted_travel_time.filled()

    # Order endpoints by distance from start
    print "Min of weighted travel time: %s, max: %s" %(np.min(weighted_travel_time),np.max(weighted_travel_time))
    print "Number of initialized endpoints is %s" % len(endx)
    Euc = []
    for i in range (0,len(endx)):
        Euc.append(np.sqrt((endx[i]-sx)**2 + (endy[i] - sy)**2 + (endz[i] - sz)**2))

    order_indici = np.argsort(Euc) # returns indices to sort elements
    Euc = np.sort(Euc)

    X = []
    Y = []
    Z = []

    for i in range(0,len(order_indici)):
        if Euc[i] > mbtr: # Check whether endpoint is sufficiently far from root voxel (min_branch_to_root)
            X.append(endx[order_indici[i]])
            Y.append(endy[order_indici[i]])
            Z.append(endz[order_indici[i]])

    print "New root is at x: %s, y: %s, z: %s" %(sx+1,sy+1,sz+1)
    print "Number of endpoints after pruning is %s" % len(X)

    sys.stdout.flush() 
    time.sleep(1)       
    # Now implement march back method to build centerline (enlarge volume)
    # The approach proceeds by updating skeleton as equal to 2
    # When branch is finished, the skeleton is solidified and set to 1
    skel = np.zeros(shape=(x_si+2,y_si+2,z_si+2), dtype = 'uint8')
    D = skel.copy() + 1.e20
    D[1:x_si+1,1:y_si+1,1:z_si+1] = weighted_travel_time
    counting = 1
    take_out = []
    number_loops = len(X)

    # Correct points for new size of volume
    start_x = sx + 1
    start_y = sy + 1
    start_z = sz + 1

    D[start_x,start_y,start_z] = 0
    skel[start_x,start_y,start_z] = 1 # initialize root

    # Begin extracting skeleton
    for ijk in range(0,number_loops):

        # Initialize endpoints and correct for larger volume   
        i = X[ijk] + 1
        j = Y[ijk] + 1
        k = Z[ijk] + 1
        
        # Check whether endpoint in neighborhood of skeleton (whisker)
        if np.all(skel[i-1:i+2,j-1:j+2,k-1:k+2])!=1: 
       
            if D[i,j,k]!=1.e20:

                done_loop = 0               
                skel[skel>0] = 1                
                
                # Check whether branch is now connected to rest of tree (stopping criteria)
                while ((i!=start_x) or (j!=start_y) or (k!=start_z)) and done_loop!=1: # can probably just do done_loop part (or) doesn't make sense (tried just done loop, always went to 1.e20 for each branch)
                #while ((i!=start_x) and (j!=start_y) and (k!=start_z)) and done_loop!=1:
                    skel[i,j,k]=2               
                    d_neighborhood = D[i-1:i+2,j-1:j+2,k-1:k+2]                    
                    
                    if np.all(skel[i-1:i+2,j-1:j+2,k-1:k+2])!=1:        
                        
                        currentMin = 1.e21 # was 1.e20
                        # Find min in neighborhood
                        for ni in range(0,3):
                            for nj in range(0,3):
                                for nk in range(0,3):
                                    if (d_neighborhood[ni,nj,nk] < currentMin) and (np.all([ni,nj,nk])!=1) and (skel[i+ni-1,j+nj-1,k+nk-1]!=2):
                                        ii = ni
                                        jj = nj
                                        kk = nk
                                        currentMin = d_neighborhood[ni,nj,nk]

                        # Update                         
                        i = i + ii - 1
                        j = j + jj - 1
                        k = k + kk - 1
                        
                        #print ijk, i,j,k, D[i,j,k]
                        #sys.stdout.flush()
                        #time.sleep(0.05)

                        if D[i,j,k] == 1.e20:
                            done_loop = 1
                            skel[skel==2] = 0 #remove branch, not marching back to root (local min in weighted_travel_time)
                    else:
                        done_loop = 1
            
            print ijk
            sys.stdout.flush() 
            time.sleep(1)                 
        else:
            take_out.append(ijk)
            

    #shift skel and start points back to correspond with original volume
    centerline_extracted = np.zeros(shape=(x_si,y_si,z_si), dtype = 'uint8')
    skel[skel==2] = 1
    centerline_extracted = skel[1:x_si+1,1:y_si+1,1:z_si+1]
    print "Number of centerline voxels is %s" %(np.sum(centerline_extracted))

    final_centerline = np.zeros(shape=(xOrig,yOrig,zOrig), dtype = 'uint8')
    final_centerline[np.min(x3):np.max(x3)+1,np.min(y3):np.max(y3)+1,np.min(z3):np.max(z3)+1] = centerline_extracted

    final_march = np.zeros(shape=(xOrig,yOrig,zOrig))
    final_march[np.min(x3):np.max(x3)+1,np.min(y3):np.max(y3)+1,np.min(z3):np.max(z3)+1] = weighted_travel_time

    return final_centerline, final_march
    

