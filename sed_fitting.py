import numpy as np
import matplotlib.pyplot as plt

magnitudes = {
    '4_14': {
        'B': {'Mag': 16.5712310415556, 'Err':  0.119056867092615},
        'r': {'Mag': 15.7421356414714, 'Err': 0.0223830455408943},
        'i': {'Mag': 15.441636338418, 'Err': 0.0434749167576849}
    },
    '4_23': {
        'r': {'Mag': 15.6703043789582, 'Err': 0.102007662572268},
        'i': {'Mag': 15.1649002694476, 'Err': 0.134509502668298}
    }
}

nu_eff = {
    'B': 4361,
    'r': 6231,
    'i': 7625
}

#First, I have to convert from apparent magnitudes to Fnu
#B filter is JC, meaning that the Magnitude-system flux zero-point is different

M_B_AB = magnitudes['4_14']['B']['Mag'] - 0.125
F_nu_B_n1 = 3631 * 10 **(-M_B_AB/2.5)
F_nu_r_n1 = 3631 * 10 **(-magnitudes['4_14']['r']['Mag']/2.5)
F_nu_i_n1 = 3631 * 10 **(-magnitudes['4_14']['i']['Mag']/2.5)
F_nu_r_n2 = 3631 * 10 **(-magnitudes['4_23']['r']['Mag']/2.5)
F_nu_i_n2 = 3631 * 10 **(-magnitudes['4_23']['i']['Mag']/2.5)

nu_B = 2.998e18 / nu_eff['B']
nu_r = 2.998e18 / nu_eff['r']
nu_i = 2.998e18 / nu_eff['i']

sigma_mag_n1 = np.array([magnitudes['4_14'][x]['Err'] for x in ['B', 'r', 'i']])
sigma_mag_n2 = np.array([magnitudes['4_23'][x]['Err'] for x in ['r', 'i']])

sigma_log_F_n1 = sigma_mag_n1 / 2.5
sigma_log_F_n2 = sigma_mag_n2 / 2.5

log_nu_n1 = np.log10([nu_B, nu_r, nu_i])
log_f_n1 = np.log10([F_nu_B_n1, F_nu_r_n1, F_nu_i_n1])
log_nu_n2 = np.log10([nu_r, nu_i])
log_f_n2 = np.log10([F_nu_r_n2, F_nu_i_n2])

weights_n1 = 1 / sigma_log_F_n1
weights_n2 = 1 / sigma_log_F_n2

(slope_n1, intercept_n1), cov_n1 = np.polyfit(log_nu_n1, log_f_n1, 1, w=weights_n1, cov='unscaled')
(slope_n2, intercept_n2), cov_n2 = np.polyfit(log_nu_n2, log_f_n2, 1, w=weights_n2, cov='unscaled')
print(f"For night 1:\nslope = {slope_n1}\nintercept = {intercept_n1}\ncovariance = {cov_n1}")
print(f"For night 2:\nslope = {slope_n2}\nintercept = {intercept_n2}\ncovariance = {cov_n2}")

def n1_fit(x): 
    y = slope_n1 * x + intercept_n1
    return y

def n2_fit(x): 
    y = slope_n2 * x + intercept_n2
    return y

x = [14.5, 14.9]
y1 = [n1_fit(x_val) for x_val in x]
y2 = [n2_fit(x_val) for x_val in x]

plt.plot(x, y1, ls='--', color='orange')
plt.plot(x, y2, ls='--', color='blue')
plt.errorbar(log_nu_n1, log_f_n1, sigma_log_F_n1, marker='o', label='4/14 Magnitudes', ls='', color='orange')
plt.errorbar(log_nu_n2, log_f_n2, sigma_log_F_n2, marker='o', label='4/23 Magnitudes', ls='', color='blue')
plt.xlabel("log(nu / Hz)")
plt.ylabel("log(F_nu / Jy)")
plt.legend()
plt.show()

#Second plot: same data and best-fit lines in linear (non-log) space
F_nu_n1 = np.array([F_nu_B_n1, F_nu_r_n1, F_nu_i_n1])
F_nu_n2 = np.array([F_nu_r_n2, F_nu_i_n2])
nu_n1 = np.array([nu_B, nu_r, nu_i])
nu_n2 = np.array([nu_r, nu_i])
sigma_F_n1 = F_nu_n1 * np.log(10) * sigma_log_F_n1
sigma_F_n2 = F_nu_n2 * np.log(10) * sigma_log_F_n2

A_n1 = 10**intercept_n1
A_n2 = 10**intercept_n2
nu_smooth = np.logspace(14.5, 14.9, 100)
F_smooth_n1 = A_n1 * nu_smooth**slope_n1
F_smooth_n2 = A_n2 * nu_smooth**slope_n2

z = 1.512
c = 2.998e18
lam_2500_obs = (1 + z) * 2500
nu_2500_obs = c / lam_2500_obs
log_nu_2500_obs = np.log10(nu_2500_obs)

log_F_2500_n1 = n1_fit(log_nu_2500_obs)
log_F_2500_n2 = n2_fit(log_nu_2500_obs)
sigma_log_F_2500_n1 = np.sqrt(log_nu_2500_obs**2 * cov_n1[0, 0] + cov_n1[1, 1] + 2 * log_nu_2500_obs * cov_n1[0, 1])
sigma_log_F_2500_n2 = np.sqrt(log_nu_2500_obs**2 * cov_n2[0, 0] + cov_n2[1, 1] + 2 * log_nu_2500_obs * cov_n2[0, 1])
F_2500_n1 = 10**log_F_2500_n1
F_2500_n2 = 10**log_F_2500_n2
sigma_F_2500_n1 = F_2500_n1 * np.log(10) * sigma_log_F_2500_n1
sigma_F_2500_n2 = F_2500_n2 * np.log(10) * sigma_log_F_2500_n2

plt.figure()
plt.errorbar(nu_2500_obs, F_2500_n1, sigma_F_2500_n1, color='orange', ls='', marker='*', ms=20, capsize=5, label='4/14 2500Å Flux')
plt.errorbar(nu_2500_obs, F_2500_n2, sigma_F_2500_n2, color='blue', ls='', marker='*', ms=20, capsize=5, label='4/23 2500Å Flux')
plt.plot(nu_smooth, F_smooth_n1, ls='--', color='orange')
plt.plot(nu_smooth, F_smooth_n2, ls='--', color='blue')
plt.errorbar(nu_n1, F_nu_n1, sigma_F_n1, marker='o', label='4/14 Magnitudes', ls='', color='orange', capsize=3)
plt.errorbar(nu_n2, F_nu_n2, sigma_F_n2, marker='o', label='4/23 Magnitudes', ls='', color='blue', capsize=3)
plt.xlabel("nu (Hz)")
plt.ylabel("F_nu (Jy)")
plt.legend()
plt.show()

#Third plot: apparent magnitude vs time for all filters and interpolated 2500AA value
from astropy.time import Time

m_2500_n1 = -2.5 * np.log10(F_2500_n1 / 3631)
m_2500_n2 = -2.5 * np.log10(F_2500_n2 / 3631)
sigma_m_2500_n1 = (2.5 / np.log(10)) * (sigma_F_2500_n1 / F_2500_n1)
sigma_m_2500_n2 = (2.5 / np.log(10)) * (sigma_F_2500_n2 / F_2500_n2)

jd_n1 = Time('2026-04-15T03:00:00', format='isot', scale='utc').jd
jd_n2 = Time('2026-04-24T02:45:00', format='isot', scale='utc').jd
jd_offset = 2460000

fig, ax = plt.subplots()
ax.errorbar(jd_n1 - jd_offset, magnitudes['4_14']['B']['Mag'], magnitudes['4_14']['B']['Err'],
            marker='o', ls='', color='steelblue', capsize=3, label='$B$')
ax.errorbar(jd_n1 - jd_offset, magnitudes['4_14']['r']['Mag'], magnitudes['4_14']['r']['Err'],
            marker='o', ls='', color='tomato', capsize=3, label='$r$')
ax.errorbar(jd_n2 - jd_offset, magnitudes['4_23']['r']['Mag'], magnitudes['4_23']['r']['Err'],
            marker='o', ls='', color='tomato', capsize=3)
ax.errorbar(jd_n1 - jd_offset, magnitudes['4_14']['i']['Mag'], magnitudes['4_14']['i']['Err'],
            marker='o', ls='', color='firebrick', capsize=3, label='$i$')
ax.errorbar(jd_n2 - jd_offset, magnitudes['4_23']['i']['Mag'], magnitudes['4_23']['i']['Err'],
            marker='o', ls='', color='firebrick', capsize=3)
ax.errorbar(jd_n1 - jd_offset, m_2500_n1, sigma_m_2500_n1,
            marker='*', ls='', ms=15, color='purple', capsize=5, label=r'$m_{2500}$ (interp.)')
ax.errorbar(jd_n2 - jd_offset, m_2500_n2, sigma_m_2500_n2,
            marker='*', ls='', ms=15, color='purple', capsize=5)
ax.invert_yaxis()
ax.set_xlabel(f"JD $-$ {jd_offset}")
ax.set_ylabel("Apparent Magnitude (AB)")
ax.legend()
plt.show()
