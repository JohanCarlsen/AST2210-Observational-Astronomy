# This program will take the spectrum file "IRAS13120_spec.txt" and fit
# a Gaussian to it.
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

_v, antennaT = np.loadtxt('IRAS13120_spec.txt', unpack=True)
v_r = _v / 1e3
S = antennaT * 42   # \pm 6

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(v_r, S, color='black' , lw=.5)

ax.set_xlabel('Radial velocity, ($v_r$) [m/s]')
ax.set_ylabel('Flux density, ($S_\\nu$) [Jy]')
plt.savefig('spectrum_plot.pdf')

def gaussian(v_r, a, b, c, d):
    '''Defining the Gaussian function to be used to fit data'''

    g = a * np.exp(-(v_r - b)**2 / (2 * c**2)) + d

    return g

b_idx = np.where(S == np.max(S))[0][0]

_a = np.max(S) - np.min(S)
_b = v_r[b_idx]
_c = .2
_d = 0

est_params = np.array([_a, _b, _c, _d])

params, covar = curve_fit(gaussian, v_r, S, est_params)
a, b, c, d = params

fitS = gaussian(v_r, a, b, c, d)
stderr = np.sqrt(np.diag(covar))

ax.plot(v_r, fitS, lw=1, color='red', label='Gauss fit')
plt.savefig('spectrum_plot_gauss.pdf')

peakX = b
peakY = np.max(fitS)

fwhm = 2 * np.sqrt(2 * np.log(2)) * c
halfMax = a / 2
fullWidth_start = b - fwhm / 2
fullWidth_stop = b + fwhm / 2
fullWidth = np.linspace(fullWidth_start, fullWidth_stop, 2)

fig, ax1 = plt.subplots(figsize=(8, 4.5))

ax1.plot(v_r, fitS, lw=1, color='red')
ax1.plot(peakX, peakY, marker='o', color='black', markersize=3)
ax1.plot(fullWidth, np.array([halfMax, halfMax]), color='black', ls='dashed', lw=1)
plt.text(peakX + .1, peakY - .05, '$S_\\nu^{peak}=$' + f'{np.max(S):.3f} $\pm$ ' + f'{stderr[0]:.3f} Jy')
plt.text(fullWidth_stop + .1, halfMax - .05, f'FWHM = {fwhm:.3f} m/s')

ax1.set_xlabel('Radial velocity, ($v_r$) [m/s]')
ax1.set_ylabel('Flux density, ($S_\\nu$) [Jy]')
ax.legend()
plt.savefig('spectrum_plot_fwhm.pdf')

plt.show()
