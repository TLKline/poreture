# Script to test centerline extraction method
# Read h61_vessel_tree file and outputs,
# testCenterline, and final_march

from poreture.extract_centerline import kline_vessel
import nibabel as nib
import numpy as np
import os

#os.chdir("..")
cur_dir = os.getcwd()
inputFILE = (cur_dir+'/h61_vessel_tree.nii.gz')
fileLOAD = nib.load(inputFILE)
fileDATA = fileLOAD.get_data()
a = fileDATA.copy()
root_location = [134,96,29]
mbtr_to_pass = 30

# Call Function
[extracted_centerline, final_march] = kline_vessel(a, root_location, min_branch_to_root = mbtr_to_pass) 
#, cluster_graph_weight = 1000, dist_map_weight = 40)

if np.any(extracted_centerline - (extracted_centerline*a)) == 1:
    print "WARNING: Centerline voxel outside original volume!!!"

# Save
output_filename = cur_dir+'/testCenterline.nii.gz'
new_image = nib.Nifti1Image(extracted_centerline,fileLOAD.get_affine())
nib.save(new_image,output_filename)

output_filename = cur_dir+'/final_march.nii.gz'
new_image = nib.Nifti1Image(final_march,fileLOAD.get_affine())
nib.save(new_image,output_filename)