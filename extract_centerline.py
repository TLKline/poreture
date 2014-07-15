#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""This module implements a number of centerline extraction methods.

The methods work on 3D binary volumes composed of
0: background
1: object to be skeletonized; have centerline extracted
The volume representations are voxel based

kline_vessel - [Kline et al. ABME 2010]

kline_pore - [Kline et al. J Porous Mat 2011]

kuba_thinning - [Kuba paper...]

This is a multi-line docstring. Paragraphs are separated with blank lines. 
Lines conform to 79-column limit. 

Module and packages names should be short, lower_case_with_underscores.

See http://www.python.org/dev/peps/pep-0008/ for more PEP-8 details and
http://wwd.ca/blog/2009/07/09/pep-8-cheatsheet/ for an up-to-date version
of this cheatsheet.
"""

# a 79-char ruler:
#234567891123456789212345678931234567894123456789512345678961234567897123456789

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

    return local_max, x, y, z

def kline_vessel(vol, startID, **kwargs):
    """This function creates a centerline from the segmented volume (vol)

    vol: 3D binary volume
    startID: index of root, as in [x,y,z] location
    
    optional kwargs:
        dist_map_weight
        cluster_graph_weight
        min_branch_length
    """

    # Imports

    # Make sure object is equal to 1
    B2 = np.array(vol.copy() > 0, dtype = 'int8')

    # Set defaults
    if kwargs is not None:
        for key, value in kwargs.iteritems():
            print "%s == %s" %(key,value)
            if key=='dist_map_weight':
                dmw = value
            if key=='cluster_graph_weight':
                cgw = value
    if 'dmw' not in locals():
        dmw = 6
    if 'cgw' not in locals():
        cgw = np.sum(vol)/20

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
    print np.min(x3),np.max(x3)
    B2 = B2[np.min(x3):np.max(x3)+1,np.min(y3):np.max(y3)+1,np.min(z3):np.max(z3)+1]

    # Setup starting index list and correct for change in volume size
    sx = startID[0] - np.min(x3) 
    sy = startID[1] - np.min(y3) 
    sz = startID[2] - np.min(z3) 

    # New volume size (bounding box)
    [x_si,y_si,z_si] = np.shape(B2)
    print x_si,y_si,z_si,sx,sy,sz

    # Perform first fast march to determine endpoints
    # works on binary speed function
    phi = B2.copy()
    phi[sx,sy,sz] = -1
    speed = np.ones(shape = (np.shape(phi)))
    mask = B2<1
    phi = np.ma.MaskedArray(phi, mask)
    binary_travel_time = skfmm.travel_time(phi, speed)

    # Fill in masked values and set to zero
    binary_travel_time = binary_travel_time.filled()
    binary_travel_time[binary_travel_time==1.e20] = 0

    # Normalize and apply cluster graph weighting
    hold_binary_travel_time = binary_travel_time.copy()
    binary_travel_time = np.round(binary_travel_time/np.max(binary_travel_time) * cgw) #by rounding group clusters together, and find maximum clusters
    #these will need to be searched then individually to find local max, one within each cluster...

    print np.max(binary_travel_time)
    # Detect local 'cluster' maximums
    [a6, endx, endy, endz] = detect_local_maxima(binary_travel_time)
    print "number of non-zero elements is %s" % (np.sum(B2))
    print "number of local maxima was %s" % (np.sum(a6))

    print np.sum(hold_binary_travel_time)
    #hold_binary_travel_time[a6==False] = 0 #remove non-max clusters
    # Not sure, but seems like maybe cluster stuff complicates it further, maybe just find local max, and deal with it at level of length of branches
    print np.sum(hold_binary_travel_time)
    [a7, endx, endy, endz] = detect_local_maxima(hold_binary_travel_time)
    print "number of local maxima was %s" % (np.sum(a7))
    return binary_travel_time, endx, endy, endz

# Import
import numpy as np
import skfmm
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology

#Test array stuff
a = np.zeros(shape=(150,150,150), dtype = 'int8')
a[5:140,5:140,5:140] = 1
a[134:136,133:141,134:136] = 1
print a.dtype

# Call Function
[travel_time, endx, endy, endz] = kline_vessel(a, [8,8,8], cluster_graph_weight = 1000, dist_map_weight = 40)

#neighborhood = morphology.generate_binary_structure(3,2)
#print neighborhood

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


   
    
