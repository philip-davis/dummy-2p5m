import struct
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

if len(sys.argv) != 2:
    print(f"{sys.argv[0]} <filename>")
    exit

#
# Binary .dat files from dummy_sim 
# -  step (4B int): timestep
# -  nx (4B int): number of data rows
# -  ny (4B int): number of data columns
# -  lb1, lb2 (8B floats): the lower bound of the data
# -  ub1, ub2 (8B floats): the upper bound of the data
# -  ps1, ps2 (8B floats): the coordinates of the plume source
# -  w1, w2 (8B floats): wind vector
# -  a (nx * ny 8B floats) particulate data
#
with open(sys.argv[1], "rb") as f:
    step, nx, ny, lb1, lb2, ub1, ub2, ps1, ps2, w1, w2 = struct.unpack("<iiidddddddd", f.read(76))
    a = np.fromfile(f)

# Switched x and y at some point...should fix this
data = a.reshape((nx, ny), order='F')


fig, ax = plt.subplots()
im = ax.imshow(data, extent = [lb1, ub1, lb2, ub2])
ax.xaxis.set_major_locator(ticker.MultipleLocator(1.00))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
plt.colorbar(im)
plt.title(f'Step {step}')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# normalization factor for the wind vector
nf = (w1**2 + w2**2)**.5

if (w1 != 0 or w2 != 0) and (ps1 != 0 or ps2 != 0):
    plt.arrow(ps1, ps2, .25*w2/h,-.25*w1/h,width=0.02, color='r')
plt.show()
