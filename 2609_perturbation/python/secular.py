import numpy as np
import sys
import matplotlib.pyplot as plt

def laplace(alpha, s, j):
    """Compute Laplace coefficient b_s^{(j)}"""
    psi = np.linspace(0, 2 * np.pi, num=100000)
    b_sj_fun = (1/np.pi) * np.cos(j * psi) / np.power(1 - 2 * alpha * np.cos(psi) + alpha * alpha, s)
    b_sj = np.trapz(b_sj_fun, psi)

    return b_sj

# Table 8.10.3 in Standish, Williams - Orbital Ephemerides of the Sun, Moon, and Planets
a_au = np.array([0.38709843, 0.72332102, 1.00000018, 1.52371243, 
                 5.20248019, 9.54149883, 19.18797948, 30.06952752])
e = np.array([0.20563661, 0.00676399, 0.01673163, 0.09336511,
              0.04853590, 0.05550825, 0.04685740, 0.00895439])
i_deg = np.array([7.00559432, 3.39777545, -0.00054346, 1.85181869,
                  1.29861416, 2.49424102, 0.77298127, 1.77005520])
varpi_deg = np.array([77.45771895, 131.76755713, 102.93005885, -23.91744784,
                      14.27495244, 92.86136063, 172.43404441, 46.68158724])
Omega_deg = np.array([48.33961819, 76.67261496, -5.11260389, 49.71320984,
                      100.29282654, 113.63998702, 73.96250215, 131.78635853])

# Table 8.8.3.
msmp = np.array([6023600, 408523.71, 328900.5614, 3098708, 
                 1047.3486, 3497.898, 22902.98, 19412.24])

# Relative mass of planets as multiples of the mass of the Sun.
M = 1 / msmp

i_rad = np.radians(i_deg)
varpi_rad = np.radians(varpi_deg)
Omega_rad = np.radians(Omega_deg)

# Define new variables.
h = e * np.sin(varpi_rad)
k = e * np.cos(varpi_rad)
p = i_rad * np.sin(Omega_rad)
q = i_rad * np.cos(Omega_rad)

num_planets = len(e)
#print(num_planets)

# The gravitational constant.
G = 6.6743e-11

# Mass of the Sun.
Msun = 1.988475e30

# Astronomical unit
au = 149597870691

# Semi-major axes in meters.
a_m = a_au * au

# Use Kepler's third law to determine the mean motion.
n = np.sqrt(G * Msun * (1.0 + M) * np.pow(a_m, -3.0))

print("n")
print(n)

# Ratio of semi-major axes.
alpha = np.zeros((num_planets, num_planets))
for ij in range(0, num_planets):
    for ik in range(0, num_planets):
        alpha[ij][ik] = a_au[ij] / a_au[ik]

# The matrices for the eigenvalue problem.
A = np.zeros((num_planets, num_planets))
B = np.zeros((num_planets, num_planets))

for ij in range(0, num_planets):
    for ik in range(0, num_planets):
        m_ratio = M[ik] / (1.0 + M[ij])
        if ij != ik:
            alpha2 = (alpha[ij][ik] ** 2)
            diag_term = 0.25 * n[ij] * m_ratio * alpha2 * laplace(alpha[ij][ik], 1.5, 1)
            A[ij][ij] = A[ij][ij] + diag_term
            B[ij][ij] = B[ij][ij] - diag_term
            A[ij][ik] = A[ij][ik] - 0.25 * m_ratio * n[ij] * alpha2 * laplace(alpha[ij][ik], 1.5, 2)
            B[ij][ik] = B[ij][ik] + 0.25 * m_ratio * n[ij] * alpha2 * laplace(alpha[ij][ik], 1.5, 1)


# Solve the eigenvalues and eigenvectors.
g, u = np.linalg.eig(A)
f, w = np.linalg.eig(B)

# Solve the scaling factors for eigenvectors and the beta/gamma angles.
M = np.block([[u, 0.0*u], [0.0*u, u]])
N = np.block([[w, 0.0*w], [0.0*w, w]])
S = np.linalg.solve(M, np.concat([h, k]))
T = np.linalg.solve(N, np.concat([p, q]))
beta = []
gamma = []

for ip in range(0, num_planets):
    value = np.degrees(np.atan2(S[ip], S[ip + num_planets]))
    value2 = np.degrees(np.atan2(T[ip], T[ip + num_planets]))
    norm = np.sqrt(S[ip] ** 2 + S[ip + num_planets] ** 2)
    norm2 = np.sqrt(T[ip] ** 2 + T[ip + num_planets] ** 2)
    if value < 0:
        value = value + 360.0
    if value2 < 0:
        value2 = value2 + 360.0
    beta.append(value)
    gamma.append(value2)
    u[:, ip] *= norm
    w[:, ip] *= norm2

np.set_printoptions(threshold=sys.maxsize,linewidth=500)
#print(np.degrees(g) * 86400 * 365.25 * 3600)
#print(np.degrees(f) * 86400 * 365.25 * 3600)

# Visualization
plt.rcParams['text.usetex'] = True
fig, ax = plt.subplots(4, 1)
fig.set_size_inches(8, 10)

# Time limits
t_min = -1000e3
t_max = 1000e3
t = np.arange(t_min, t_max, 100)
yeartosec = 365.25 * 86400
tsec = yeartosec * t

ht = np.zeros((num_planets, len(t)))
kt = np.zeros((num_planets, len(t)))
pt = np.zeros((num_planets, len(t)))
qt = np.zeros((num_planets, len(t)))

for planet in range(0, num_planets):
    for ip in range(0, num_planets):
        ht[planet, :] = ht[planet, :] + u[planet, ip] * np.sin(g[ip] * tsec + np.radians(beta[ip]))
        kt[planet, :] = kt[planet, :] + u[planet, ip] * np.cos(g[ip] * tsec + np.radians(beta[ip]))
        pt[planet, :] = pt[planet, :] + w[planet, ip] * np.sin(f[ip] * tsec + np.radians(gamma[ip]))
        qt[planet, :] = qt[planet, :] + w[planet, ip] * np.cos(f[ip] * tsec + np.radians(gamma[ip]))

ecc = np.sqrt(ht * ht + kt * kt)
varpi = np.degrees(np.atan2(ht, kt))
I = np.degrees(np.sqrt(pt * pt + qt * qt))
Omega = np.degrees(np.atan2(pt, qt))

styles = ['k', 'r:', 'b', 'g:', 'k', 'r:', 'b', 'g:']
names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
line_widths = [0.5, 0.5, 1, 1, 0.5, 0.5, 1, 1]

for planet in range(4, 8):
    ax[0].plot(t / 1000.0, ecc[planet, :], styles[planet], label=names[planet], linewidth=line_widths[planet])
    ax[1].plot(t / 1000.0, varpi[planet, :], styles[planet], label=names[planet], linewidth=line_widths[planet])
    ax[2].plot(t / 1000.0, I[planet, :], styles[planet], label=names[planet], linewidth=line_widths[planet])
    ax[3].plot(t / 1000.0, Omega[planet, :], styles[planet], label=names[planet], linewidth=line_widths[planet])

ax[0].set(xlabel='', ylabel="e")
ax[0].set_xlim(t_min/1000.0, t_max/1000.0)
ax[0].grid()
ax[0].legend(loc="lower right")
ax[1].set(xlabel='', ylabel=r"$\varpi$ (deg)")
ax[1].set_xlim(t_min/1000.0, t_max/1000.0)
ax[1].grid()
ax[2].set(xlabel='', ylabel="i (deg)")
ax[2].set_xlim(t_min/1000.0, t_max/1000.0)
ax[2].grid()
ax[3].set(xlabel='Time (millenia)', ylabel=r"$\Omega$ (deg)")
ax[3].set_xlim(t_min/1000.0, t_max/1000.0)
ax[3].grid()

fig.savefig('fig/major_planets.eps', format='eps')
