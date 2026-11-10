import numpy as np
import matplotlib.pyplot as plt

t_min = -1e5
t_max = 1e5
t = np.arange(t_min, t_max, 100)
beta1_deg = -146.89233296
beta2_deg = -53.35651935
gamma1_deg = 105.73951730
gamma2_deg = 126.82475444
beta1_rad = np.radians(beta1_deg)
beta2_rad = np.radians(beta2_deg)
gamma1_rad = np.radians(gamma1_deg)
gamma2_rad = np.radians(gamma2_deg)

g_1_deg = 0.0009634346412059512
g_2_deg = 0.006099080354401151
f_1_deg = 0
f_2_deg = -0.007062514995607103
g_1_rad = np.radians(g_1_deg)
g_2_rad = np.radians(g_2_deg)
f_1_rad = np.radians(f_1_deg)
f_2_rad = np.radians(f_2_deg)

u_11 = -0.04388212 
u_21 = -0.03543747
u_12 = 0.01557876 
u_22 = -0.04758095
w_11 = 0.02856889 
w_21 = 0.02856889
w_12 = -0.00629765
w_22 = 0.01553295

h_1 = u_11 * np.sin(g_1_rad * t + beta1_rad) + u_12 * np.sin(g_2_rad * t + beta2_rad)
h_2 = u_21 * np.sin(g_1_rad * t + beta1_rad) + u_22 * np.sin(g_2_rad * t + beta2_rad)
k_1 = u_11 * np.cos(g_1_rad * t + beta1_rad) + u_12 * np.cos(g_2_rad * t + beta2_rad)
k_2 = u_21 * np.cos(g_1_rad * t + beta1_rad) + u_22 * np.cos(g_2_rad * t + beta2_rad)

p_1 = w_11 * np.sin(f_1_rad * t + gamma1_rad) + w_12 * np.sin(f_2_rad * t + gamma2_rad)
p_2 = w_21 * np.sin(f_1_rad * t + gamma1_rad) + w_22 * np.sin(f_2_rad * t + gamma2_rad)
q_1 = w_11 * np.cos(f_1_rad * t + gamma1_rad) + w_12 * np.cos(f_2_rad * t + gamma2_rad)
q_2 = w_21 * np.cos(f_1_rad * t + gamma1_rad) + w_22 * np.cos(f_2_rad * t + gamma2_rad)

ecc_1 = np.sqrt(h_1 * h_1 + k_1 * k_1)
ecc_2 = np.sqrt(h_2 * h_2 + k_2 * k_2)
varpi_1 = np.degrees(np.atan2(h_1, k_1))
varpi_2 = np.degrees(np.atan2(h_2, k_2))
I_1 = np.degrees(np.sqrt(p_1 * p_1 + q_1 * q_1))
I_2 = np.degrees(np.sqrt(p_2 * p_2 + q_2 * q_2))
Omega_1 = np.degrees(np.atan2(p_1, q_1))
Omega_2 = np.degrees(np.atan2(p_2, q_2))

fig, ax = plt.subplots(2, 2)
ax[0][0].plot(t/1000.0, ecc_1, 'k', label="Jupiter")
ax[0][0].plot(t/1000.0, ecc_2, 'k:', label="Saturn")
ax[0][0].grid()
ax[0][0].set_xlim(t_min/1000.0, t_max/1000.0)
ax[0][0].set(xlabel='Time (millenia)', ylabel="e")
ax[0][0].legend()

ax[0][1].plot(t/1000.0, varpi_1, 'k', label="Jupiter")
ax[0][1].plot(t/1000.0, varpi_2, 'k:', label="Saturn")
ax[0][1].grid()
ax[0][1].set_xlim(t_min/1000.0, t_max/1000.0)
ax[0][1].set(xlabel='Time (millenia)', ylabel="ϖ (deg)")
#ax[0][1].legend()

ax[1][0].plot(t/1000.0, I_1, 'k', label="Jupiter")
ax[1][0].plot(t/1000.0, I_2, 'k:', label="Saturn")
ax[1][0].grid()
ax[1][0].set_xlim(t_min/1000.0, t_max/1000.0)
ax[1][0].set(xlabel='Time (millenia)', ylabel="I (deg)")
#ax[1][0].legend()

ax[1][1].plot(t/1000.0, Omega_1, 'k', label="Jupiter")
ax[1][1].plot(t/1000.0, Omega_2, 'k:', label="Saturn")
ax[1][1].grid()
ax[1][1].set_xlim(t_min/1000.0, t_max/1000.0)
ax[1][1].set(xlabel='Time (millenia)', ylabel="Ω (deg)")
#ax[1][1].legend()

fig.set_size_inches(12, 7)

plt.show()

fig.savefig('fig/twoplanet_elements.eps', format='eps')