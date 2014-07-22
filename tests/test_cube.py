# Test on 3D cube
from poreture.extract_centerline import kline_vessel
import numpy as np

# Create cube
a = np.zeros(shape=(150,150,150), dtype = 'int8')
a[5:140,5:140,5:140] = 1
a[134:136,133:145,134:136] = 1
root_location = [8,8,8]
mbtr_to_pass = 20

# Call Function
[extracted_centerline, final_march] = kline_vessel(a, root_location, min_branch_to_root = mbtr_to_pass) 