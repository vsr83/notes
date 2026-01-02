import numpy as np
import matplotlib.pyplot as plt

from rgb import cie_xyz_integrate, cie_xyz_delta, cie_xyz_sgrb
from blackbody import blackbody_spectrum_wavel

SRGB = []

for T in range(800, 12000):
    wavel_nm, spectrum = blackbody_spectrum_wavel(T, 360, 830, 1)
    X, Y, Z = cie_xyz_integrate(wavel_nm, spectrum)
    M = max(X, Y, Z)
    R, G, B = cie_xyz_sgrb(X/M, Y/M, Z/M)

    SRGB.append([R, G, B])

    #print("%f %f %f" % (X, Y, Z))
    #print("%f %f %f" % (R, G, B))

SRGB = np.array([SRGB])
fig, ax = plt.subplots()
ax.imshow(SRGB, extent=[800, 12000, 0, 1000])

ax.set_yticks([])
#ax.set_xticks([1000, 3000, 5000, 7000, 9000, 11000, 12000])
ax.set_xlabel("Effective Temperature (K)")

fig.set_size_inches(7, 2)
fig.savefig('../fig/effective_temperature.eps', format='eps')

fig2, ax2 = plt.subplots()
wavel_nm, spectrum1 = blackbody_spectrum_wavel(4000, 250, 2500, 1)
wavel_nm, spectrum2 = blackbody_spectrum_wavel(5000, 250, 2500, 1)
wavel_nm, spectrum3 = blackbody_spectrum_wavel(6000, 250, 2500, 1)
ax2.plot(wavel_nm, spectrum1, label="4000 K")
ax2.plot(wavel_nm, spectrum2, label="5000 K")
ax2.plot(wavel_nm, spectrum3, label="6000 K")
ax2.set_ylabel("Spectral Radiance (W/m^2/nm)")
ax2.set_xlabel("Wavelength (nm)")
ax2.grid()
ax2.legend()
fig2.savefig('../fig/blackbody_example.eps', format='eps')


fig3, ax3 = plt.subplots()
SRGB2 = []
for beta in np.arange(-0.99, 1.00, 0.01):
    print(beta)
    SRGB = []
    for T in np.arange(800, 12000, 10):
        Teff = T * np.sqrt((1 + beta) / (1 - beta))

        wavel_nm, spectrum = blackbody_spectrum_wavel(Teff, 360, 830, 1)
        X, Y, Z = cie_xyz_integrate(wavel_nm, spectrum)
        M = max(X, Y, Z)
        R, G, B = cie_xyz_sgrb(X/M, Y/M, Z/M)

        SRGB.append([R, G, B])
    SRGB2.append(SRGB)


SRGB2 = np.flip(np.array(SRGB2), axis=0)
ax3.imshow(SRGB2, extent=[800, 12000, 0, 10000])

ax3.set_yticks([0, 5000, 10000], ['-1.0', '0.0', '1.0'])
#ax.set_xticks([1000, 3000, 5000, 7000, 9000, 11000, 12000])
ax3.set_xlabel("Effective Temperature (K)")
ax3.set_ylabel("Observer Velocity (β)")
#ax3.grid(color='k')

fig4, ax4 = plt.subplots()
SRGB3 = []
T = 5778
for beta in np.arange(-0.99, 1.00, 0.01):
    print(beta)
    SRGB = []
    gamma = 1/np.sqrt(1 - beta*beta)
    for theta in np.arange(-90, 90, 0.1):
        theta_rad = np.radians(theta)
        Teff = T / gamma / (1 - beta * np.cos(theta_rad))

        wavel_nm, spectrum = blackbody_spectrum_wavel(Teff, 360, 830, 1)
        X, Y, Z = cie_xyz_integrate(wavel_nm, spectrum)
        M = max(X, Y, Z)
        R, G, B = cie_xyz_sgrb(X/M, Y/M, Z/M)

        SRGB.append([R, G, B])
    SRGB3.append(SRGB)


SRGB3 = np.flip(np.array(SRGB3), axis=0)
ax4.imshow(SRGB3, extent=[0, 10000, 0, 10000])

beta = np.arange(0, 1, 0.01)
gamma = 1/np.sqrt(1 - beta*beta)
theta = np.degrees(np.acos(1/beta * (1 - 1/gamma)))
ax4.plot(5000 * (1 + theta/90), 5000 * (1 + beta), 'k:')
ax4.plot(5000 * (1 - theta/90), 5000 * (1 + beta), 'k:')
ax4.plot(np.array([0, 10000]), np.array([5000, 5000]), 'k:')

ax4.set_yticks([0, 5000, 10000], ['-1', '0.0', '1'])
ax4.set_xticks([0, 2500, 5000, 7500, 10000], ['-90', '-45', '0', '45', '90'])
ax4.set_xlabel("Angle to Target (θ)")
ax4.set_ylabel("Observer Velocity (β)")
ax4.set_title("$\mathregular{T_{eff} = 5778} K$")
#ax3.grid(color='k')
plt.show()