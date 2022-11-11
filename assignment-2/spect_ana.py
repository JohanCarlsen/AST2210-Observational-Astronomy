# This program will take the spectrum file "IRAS13120_spec.txt" and fit
# a Gaussian to it.
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

v_rClass, antennaT = np.loadtxt('IRAS13120_spec.txt', unpack=True)
v_rGauss = np.linspace(v_rClass[0], v_rClass[-1], 1001)
S = antennaT * 42   # \pm 6

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(v_rClass, S, color='black' , lw=.5)

ax.set_xlabel('Radial velocity, ($v_r$) [km/s]')
ax.set_ylabel('Flux density, ($S_\\nu$) [Jy]')
plt.savefig('spectrum_plot.pdf')

def gaussian(v_r, a, b, c, d):
    '''Defining the Gaussian function to be used to fit data'''

    g = a * np.exp(-(v_r - b)**2 / (2 * c**2)) + d

    return g

b_idx = np.where(S == np.max(S))[0][0]

_a = np.max(S) - np.min(S)
_b = v_rClass[b_idx]
_c = 200
_d = 0

v_RMSClass = 3.074e-3 * 42  # pm 6

est_params = np.array([_a, _b, _c, _d])
# print(est_params)
'''
[  3.17898   5.533   200.        0.     ]
'''

params, covar = curve_fit(gaussian, v_rClass, S, est_params)
a, b, c, d = params
# print(params)
'''
[ 2.85571042e+00  8.17360151e+01  1.36828958e+02 -9.30479244e-03]
'''

fitS = gaussian(v_rGauss, a, b, c, d)
stderr = np.sqrt(np.diag(covar))
# print(stderr)
'''
[0.05750111 3.15418348 3.23577549 0.01066059]
'''

ax.plot(v_rGauss, fitS, lw=1, color='red', label='Gauss fit')
ax.legend()
plt.savefig('spectrum_plot_gauss.pdf')

peakX = b
peakY = a - abs(d)
# print(f'S/N = {peakY / v_RMSClass:.3f}')
'''
S/N = 22.047
'''

# Error on the total line flux
dFa = c * np.sqrt(2 * np.pi)
dA = stderr[0]
dFc = (a-abs(d)) * np.sqrt(2 * np.pi)
dC = stderr[2]

dF = np.sqrt(dFa**2 * dA**2 + dFc**2 + dC**2)
# print(dF)
'''
21.220789498985443
'''

fwhm = 2 * np.sqrt(2 * np.log(2)) * c
halfMax = a / 2
fullWidth_start = b - fwhm / 2
fullWidth_stop = b + fwhm / 2
fullWidth = np.linspace(fullWidth_start, fullWidth_stop, 2)

# Error in the FWHM
dFWHM = 2 * np.sqrt(2 * np.log(2)) * dC
# print(dFWHM)
'''
7.619668975054973
'''

fig, ax1 = plt.subplots(figsize=(8, 4.5))

ax1.plot(v_rGauss, fitS, lw=1, color='red')
ax1.plot(peakX, peakY, marker='o', color='black', markersize=3)
ax1.plot(fullWidth, np.array([halfMax, halfMax]), color='black', ls='dashed', lw=1)
plt.text(peakX + 100, peakY - .05, '$S_\\nu^{peak}=$' + f'{a-abs(d):.3f} $\pm$ ' + f'{stderr[0]:.3f} Jy')
plt.text(fullWidth_stop + 100, halfMax - .05, f'FWHM = {fwhm:.3f} $\pm$ ' + f'{dFWHM:.3f} km/s')

ax1.set_xlabel('Radial velocity, ($v_r$) [km/s]')
ax1.set_ylabel('Flux density, ($S_\\nu$) [Jy]')
ax.legend()
plt.savefig('spectrum_plot_fwhm.pdf')

plt.show()
