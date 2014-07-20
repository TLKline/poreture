# -*- coding: utf-8 -*-

"""
This module implements centerline extraction and skeletonization methods.

The methods work on 3D binary volumes composed of
    0: background
    1: object to be skeletonized/centerline extracted

The work was inspired by the need for analysis of voxel-based geometries

Literature references for implemented methods
    kline_vessel - [Kline et al. ABME 2010]
    kline_pore - [Kline et al. J Porous Mat 2011]
    kuba_thinning - [Kuba paper...]
"""

# a 79-char ruler:
#234567891123456789212345678931234567894123456789512345678961234567897123456789

# Imports
import os
import sys
import time
import nibabel as nib
import numpy as np
import skfmm
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

def find_3D_object_voxel_list(vol):
    
    B2 = np.array(vol.copy() > 0)

def find_terminal_end_points(vol):

    something = 1

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

    # Find 3D coordinates of volume (in order to limit size of volume)
    # There has to be a much faster, easier, prettier way!!!
    x3 = []
    y3 = []
    z3 = []
    for i in range(0,xOrig):
        for j in range(0,yOrig):
            for k in range(0,zOrig):
                if B2[i,j,k] > 0:
                    x3.append(i)
                    y3.append(j)
                    z3.append(k)

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
    

test_case = 2

#Test array for development
if test_case == 1:
    a = np.zeros(shape=(150,150,150), dtype = 'int8')
    a[5:140,5:140,5:140] = 1
    a[134:136,133:145,134:136] = 1
    root_location = [8,8,8]
    mbtr_to_pass = 20

if test_case == 2:
    cur_dir = os.getcwd()
    inputFILE = (cur_dir+'/DATA/'+'h61cu_resized_t400.nii.gz')
    fileLOAD = nib.load(inputFILE)
    fileDATA = fileLOAD.get_data()
    a = fileDATA.copy()
    root_location = [134,96,29]
    mbtr_to_pass = 30
    print np.shape(a)
    print fileLOAD.get_affine()
    print a[96,54,30]
    print np.min(a),np.max(a)


# Call Function
[extracted_centerline, final_march] = kline_vessel(a, root_location, min_branch_to_root = mbtr_to_pass) 
#, cluster_graph_weight = 1000, dist_map_weight = 40)

if np.any(extracted_centerline - (extracted_centerline*a)) == 1:
    print "WARNING: Centerline voxel outside original volume!!!"

# Save
output_filename = cur_dir+'/DATA/testCenterline4.nii.gz'
new_image = nib.Nifti1Image(extracted_centerline,fileLOAD.get_affine())
nib.save(new_image,output_filename)

output_filename = cur_dir+'/DATA/final_march4.nii.gz'
new_image = nib.Nifti1Image(final_march,fileLOAD.get_affine())
nib.save(new_image,output_filename)

"""
phi = np.ones((3,3,3)) 
phi[1,1,1] = -1
phi[0,1,1] = 0
phi[1,2,1] = 0
speed = np.ones((3,3,3)) * 200 #faster speed equals faster travel time (i.e., lower return value from skfmm.travel_time)
print skfmm.travel_time(phi, speed)

    # Perform binary fmm to determine terminal endpoints

%fmm for binary (no weight function)
tic;[D] = fast_marching(double(B2),[start_x;start_y;start_z]);toc; %was t7 before
%save -v7.3 D D

%begin determining endpoints
b = D;
b(D==Inf)=0;
[xxx,yyy,zzz]=find3d_max_point(b);
b = b./b(xxx(1),yyy(1),zzz(1));
b = round(b.*cluster_ref);
a6 = imregionalmax(b);
clear b xxx yyy zzz
[L,N]=bwlabeln(a6);
clear a6
N
%put this stuff back in
a4 = bwdist_old(~B2);
%save -v7.3 a4 a4

%take this out
%load a4 a4

%new_data = a4.*a4; %was just a4 times some constant
done_bw = 1
clear B2

[X_m,Y_m,Z_m]=find3d_max_point(a4);
max_val = a4(X_m(1),Y_m(1),Z_m(1));
a4=a4./max_val;


tic
stats = regionprops(L,D,'PixelList','PixelValues');

%save stats stats

s = size(stats);s = s(1,1)
for i = 1:s; 
    [mmm2,iii2] = max(stats(i,1).PixelValues);
    Y(i) = round(stats(i,1).PixelList(iii2,1));
    X(i) = round(stats(i,1).PixelList(iii2,2));
    Z(i) = round(stats(i,1).PixelList(iii2,3));
end
toc


while 0
tic
   for i = 1:N
        %ss = (L==i);
        f = find(L==i);
        [mm,ii] = max(D(f));
        [x3,y3,z3] = ind2sub(size(D),f(ii));
        X(i)=x3(1);
        Y(i)=y3(1);
        Z(i)=z3(1);
   end
end


%save endPoints22 X Y Z
   
   done_finding_endpoint = 1
    size(X)

%clear memory
clearvars -except a4 dist_map_weight start_x start_y start_z X Y Z header output_filename x_si y_si z_si x3 y3 z3 xOrig yOrig zOrig

%2nd fmm with weighting
tic;[D3] = fast_marching(a4.^dist_map_weight,[start_x;start_y;start_z]);toc;
end_marching = 1

%order endpoints by distance from start
for i = 1:numel(X)
    Euc(i) = sqrt((X(i)-start_x)^2 + (Y(i) - start_y)^2 + (Z(i) - start_z)^2);
end
[s,swapping] = sort(Euc);
swapping = fliplr(swapping);
X = X(swapping);
Y = Y(swapping);
Z = Z(swapping);

nn22 = sum(Euc<60) %was 30
if nn22~=0
    X(end-nn22+1:end) = [];
    Y(end-nn22+1:end) = [];
    Z(end-nn22+1:end) = [];
end

%free up more memory
clear new2 B2 
clear B2 a4 a6 a7 L b
clear a4 a6 a7 ss X_m2 Y_m2 Z_m2 X_m Y_m Z_m D

number_loops = numel(X)

%compensate for larger volume
start_x = start_x + 1;
start_y = start_y + 1;
start_z = start_z + 1;
X = X + 1;
Y = Y + 1;
Z = Z + 1;

skel = zeros(x_si+2,y_si+2,z_si+2);
skel(start_x,start_y,start_z)=1;
skel = uint8(skel);

%[xx22,yy22,zz22] = size(skel);

D_hold = Inf(size(skel));
D_hold(2:end-1,2:end-1,2:end-1) = D3;
D3 = D_hold;
clear D_hold

counting = 1;
take_out = [];

print #number of endpoints found

%begin extracting skeleton
for ijk = 1:number_loops
    
    i = X(ijk);
    j = Y(ijk);
    k = Z(ijk);
    
    
    if skel(i-1:i+1,j-1:j+1,k-1:k+1)~=1   %was 2's
   
    if D3(i,j,k)~=Inf
    done_loop = 0;
   
    skel(skel>0) = 1;
    
    
    while ((i~=start_x)||(j~=start_y)||(k~=start_z))&&done_loop~=1
   
        skel(i,j,k)=2;
       
    
        d_neighborhood = D3(i-1:i+1,j-1:j+1,k-1:k+1);
        
        
        
        
        if skel(i-1:i+1,j-1:j+1,k-1:k+1)~=1        
            
        f = find(d_neighborhood == min(d_neighborhood(:)));
        [ii,jj,kk] = ind2sub(size(d_neighborhood),f);    
        
        i = i + ii(1) - 2;
        j = j + jj(1) - 2;
        k = k + kk(1) - 2;
       
        
        else
            done_loop = 1;
                      
        end

       
       
    end
    end
    else
        take_out(counting) = ijk;
        counting = counting + 1
    end
    ijk
   
end
toc

%shift skel and start points back to correspond with original volume
skel(skel==2) = 1;
skel2 = skel(2:end-1,2:end-1,2:end-1);
clear skel
start_x = start_x - 1;start_y = start_y - 1;start_z = start_z - 1;

%make sure one connected object
CC = bwconncomp2(skel2);
largest = 0; 
for i = 1:CC.NumObjects; 
    if numel(CC.PixelIdxList{1,i})>largest; 
        indexing = i; 
        largest = numel(CC.PixelIdxList{1,i});
    end; 
end;
CC2.PixelIdxList{1,1} = CC.PixelIdxList{1,indexing};
CC2.Connectivity = 26;
CC2.NumObjects = 1;
CC2.ImageSize = CC.ImageSize;
skel2 = labelmatrix(CC2)>0;

%remove any unnecesarry points
skel2 = kuba_newThin_t3(skel2);

clearvars -except skel2 header output_filename x_si y_si z_si x3 y3 z3 xOrig yOrig zOrig start_x start_y start_z
out = zeros(xOrig,yOrig,zOrig,'uint8');
out(min(x3):max(x3),min(y3):max(y3),min(z3):max(z3)) = skel2;
clear skel2 

start_x = start_x + min(x3) - 1
start_y = start_y + min(y3) - 1
start_z = start_z + min(z3) - 1
%save file
flag2 = writeavw_original8(uint8(out),header,output_filename);

"""


   
    
