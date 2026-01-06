import numpy as np
import matplotlib.pyplot as plt 

from rgb import cie_xyz_integrate, cie_xyz_delta, cie_xyz_sgrb

SRGB = []
SRGB2 = []

for wavel in np.arange(380, 700):
    X, Y, Z = cie_xyz_delta(wavel)
    M = max(X, Y, Z)
    R, G, B = cie_xyz_sgrb(X/M, Y/M, Z/M)
    R2, G2, B2 = cie_xyz_sgrb(X, Y, Z)
    SRGB.append([R, G, B])
    SRGB2.append([R2, G2, B2 ])

fig = plt.figure()
ax0 = fig.add_subplot(211)
ax1 = fig.add_subplot(212)
SRGB = np.array([SRGB])
SRGB2 = np.array([SRGB2])

ax0.imshow(SRGB, extent=[380, 700, 0, 50])
ax0.set_yticks([])
ax0.set_xticks([380, 400, 450, 500, 550, 600, 650, 700])
ax0.set_xlabel("Wavelength (nm)")
ax0.set_title("Normalized")

ax1.imshow(SRGB2, extent=[380, 700, 0, 50])
ax1.set_yticks([])
ax1.set_xticks([380, 400, 450, 500, 550, 600, 650, 700])
ax1.set_xlabel("Wavelength (nm)")
ax1.set_title("Constant-Amplitude")

fig.set_size_inches(8, 4)
fig.savefig('../fig/rgb.eps', format='eps')

plt.show()