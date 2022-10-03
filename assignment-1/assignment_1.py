# Code for assignment 1
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import seaborn as sns
import scipy.constants as const
from numba import njit

sns.color_palette('bright')
sns.set_theme()
sns.set_style('darkgrid')

plt.rcParams['image.origin'] = 'lower'
plt.rcParams['image.cmap'] = 'viridis'

'''Loading data'''
idata = np.load('idata_square.npy')     # Spectrum data (y,x,lambda)
spect_pos = np.load('spect_pos.npy')    # Wavelength in [Å]

'''Points to observe'''
A = np.array([49, 197])
B = np.array([238, 443])
C = np.array([397, 213])
D = np.array([466, 52])

obs_points = np.array([A, B, C, D])
obs_names = np.array(['A', 'B', 'C', 'D'])
colors = np.array(['black', 'red', 'orange', 'blue'])

'''Sub field of view'''
x, y = (525, 325)
h = 100
w = 150

rect = Rectangle((x-1, y-1), w, h, linewidth=2, edgecolor='yellow',
facecolor='none')

x1, y1 = rect.get_xy()
x2 = x1 + rect.get_width()
y2 = y1 + rect.get_height()

slice_x = slice(x1, x2)
slice_y = slice(y1, y2)

idata_cut = idata[slice_y, slice_x, :]


def plot_intensity(wav_idx, include_points=None, add_patch=None, save=None):
    '''
    Plot the intensity of given wavelength
    from the entire range
    '''

    wav_value = spect_pos[wav_idx]
    intensity_data = idata[:, :, wav_idx]

    fig, ax = plt.subplots()
    ax.grid(False)

    im = ax.imshow(intensity_data)
    cbar = fig.colorbar(im)

    if include_points == True:

        for i in range(len(obs_points)):

            x, y = obs_points[i]
            circ = Circle((x-1, y-1), radius=10, facecolor='none',
            edgecolor='yellow', linewidth=1)

            ax.add_patch(circ)
            plt.text(x-1, y+15, obs_names[i], horizontalalignment='center',
            fontsize='large', color='yellow')

    cbar.ax.set_ylabel(r'Intensity')
    ax.set_title(f'Intensity for $\lambda_i$ = {wav_value} Å')
    ax.set_xlabel('x [idx]')
    ax.set_ylabel('y [idx]')

    fig.tight_layout()

    if add_patch == True:

        ax.add_patch(rect)

    if save == True:

        plt.savefig('intensity.pdf')
        plt.show()


def plot_intensity_sub(wav_idx, include_points=None, save=None):
    '''Plot the intensity of the sub field of view'''

    wav_value = spect_pos[wav_idx]
    intensity_data = idata_cut[:, :, wav_idx]

    fig, ax = plt.subplots()
    ax.grid(False)

    im = ax.imshow(intensity_data)
    cbar = fig.colorbar(im)

    if include_points == True:

        for i in range(len(obs_points)):

            x, y = obs_points[i]
            circ = Circle((x-1, y-1), radius=10, facecolor='none',
            edgecolor='yellow', linewidth=1)

            ax.add_patch(circ)
            plt.text(x-1, y+15, obs_names[i], horizontalalignment='center',
            fontsize='large', color='yellow')

    cbar.ax.set_ylabel(r'Intensity')
    ax.set_title(f'Intensity for $\lambda_i$ = {wav_value} Å of the sub FoV')
    ax.set_xlabel('x [idx]')
    ax.set_ylabel('y [idx]')

    fig.tight_layout()

    if save == True:

        plt.savefig('intensity_sub_fov.pdf')
        plt.show()


def plot_spectrum(point, obs_name, add_gaussian=None, add_average=None,
save=None):
    '''
    Take a point (x,y) in pixels and plot the
    intensity as a function of wavelength
    '''

    x, y = point

    wavelength_spectrum = idata[y-1, x-1, :]

    fig, ax = plt.subplots(figsize=(7, 5))

    if add_average == True:

        avg1 = np.average(idata, axis=0)
        avg2 = np.average(avg1, axis=0)

        ax.plot(spect_pos, avg2, ls='dashed', lw=1, color='royalblue',
        label='Average', marker='x')


    ax.plot(spect_pos, wavelength_spectrum, ls='dashed', lw=1, color='red',
    marker='x', label=f'Intensity')

    if add_gaussian == True:

        parameters, x_fit, g_fit = fit_gaussian(point)
        ax.plot(x_fit, g_fit, color='k', lw=1, label='Gauss fit')

    ax.set_title(f'Spectrum at {obs_name} = ({x},{y}) px')
    ax.set_xlabel(r'$\lambda$ in [Å]')
    ax.set_ylabel('Intensity')
    ax.legend()

    fig.tight_layout()


    if save == True:

        if add_gaussian == True:

            plt.savefig(f'spectre_at_{obs_name}_with_gauss_fit.pdf')
            plt.show()

        else:

            plt.savefig(f'spectre_at_{obs_name}.pdf')
            plt.show()


def gaussian(wavelength, a, b, c, d):
    '''Defining the Gaussian function to be used to fit data'''

    g = a * np.exp(-(wavelength - b)**2 / (2 * c)**2) + d

    return g


def fit_gaussian(point=None, average=None):
    '''
    The function takes x and y coords from point argument and
    retrieves the corresponding spectrum from the dataset.
    The parameters a, b, c and d are calculated and used in
    scipy.optimize.curve_fit to get new parameters that will be used in the
    function gaussian.
    '''

    if average == True:

        avg1 = np.average(idata, axis=0)
        avg2 = np.average(avg1, axis=0)
        wav_spec = avg2

    else:

        x, y = point
        wav_spec = idata[y-1, x-1, :]

    b_index = np.where(wav_spec == min(wav_spec))[0][0]

    _a = np.min(wav_spec) - np.max(wav_spec)
    _b = spect_pos[b_index]
    _c = .05
    _d = np.max(wav_spec)

    est_params = np.array([_a, _b, _c, _d])

    parameters, covariance = curve_fit(gaussian, spect_pos, wav_spec, est_params)
    a, b, c, d = parameters

    xaxis_fine = np.linspace(spect_pos[0], spect_pos[-1], 101)

    g = gaussian(xaxis_fine, a, b, c, d)

    return parameters, xaxis_fine, g


def estimate_central_wavelength(plot=True):
    '''Estimate the central wavelength of Fe I'''

    parameters, x_fit, g_fit = fit_gaussian(average=True)

    central_wavelength = parameters[1]

    if plot == False:

        return central_wavelength

    fig, ax = plt.subplots()

    ax.plot(x_fit, g_fit, lw=1, color='red', label='Gaussian fit')
    ax.plot(parameters[1], np.min(g_fit), 'x', color='black')

    ax.set_title('Averaged wavelength')
    ax.set_xlabel(r'$\lambda$ in [Å]')
    ax.set_ylabel('Intensity')
    ax.set_ylim(225, 450)

    plt.text(central_wavelength, np.min(g_fit)-15,
    f'$({central_wavelength:.4f},{np.min(g_fit):.4f})$',
    horizontalalignment='center')

    ax.legend()
    fig.tight_layout()

    plt.savefig('central_wavelength.pdf')


def doppler_velocity(lambda_obs_in_Å):
    '''Calculate the doppler velocity for a given wavelength'''

    Å = const.angstrom  # Angstrom
    c = const.c         # Speed of light in m/s

    lambda_em = estimate_central_wavelength(plot=False) * Å
    lambda_obs = lambda_obs_in_Å * Å

    radial_velocity = c * (lambda_obs - lambda_em) / lambda_em

    return radial_velocity


def doppler_velocity_point(point):
    '''Calculate the doppler velocity at a given point (x,y)'''

    params, x, g = fit_gaussian(point)
    lambda_obs = params[1]

    vel_r = doppler_velocity(lambda_obs)

    return vel_r


def doppler_map(title, field_of_view=None, save=None):
    '''For anym point (x,y) in the fov, make a gaussian, calculate the central
    wavelength and compare to the average. Then make a plot of the results'''

    if field_of_view == 'sub':

        M = len(idata_cut[:, 0, 0])
        N = len(idata_cut[0, :, 0])

        start_x = 525
        start_y = 325
        stop_x = start_x + 150
        stop_y = start_y + 100

    else:

        M = len(idata[:, 0, 0])
        N = len(idata[0, :, 0])

        start_x = 0
        start_y = 0
        stop_x = 750
        stop_y = 550

    velocities = np.zeros((M, N))

    for x in range(start_x, stop_x):

        for y in range(start_y, stop_y):

            point = (x, y)

            vel_r = doppler_velocity_point(point)

            velocities[y - start_y, x - start_x] = vel_r

    fig, ax = plt.subplots()
    ax.grid(False)

    im = ax.imshow(velocities, cmap='seismic')
    cbar = fig.colorbar(im)

    cbar.ax.set_ylabel(r'Velocity')
    ax.set_title(f'Doppler map for {title}')
    ax.set_xlabel('x [idx]')
    ax.set_ylabel('y [idx]')

    fig.tight_layout()

    if save == True:

        plt.savefig(f'doppler_map_{title}.pdf')


'''Plotting spectra for points A, B, C and D'''
# for i in range(4):
#
#     plot_spectrum(obs_points[i], obs_names[i], add_average=True, save=True)

'''Plotting the spectre at the points A, B, C and D together
with the fittet Gauss curve'''
# for i in range(4):
#
#     plot_spectrum(obs_points[i], obs_names[i], add_gaussian=True, save=True)

'''Estimating the central wavelength of Fe I'''
# cent_wav = estimate_central_wavelength()

'''Estimate doppler velocity for points A, B, C and D'''
# print('\n\tPoint\tDoppler velocity\n')
#
# for i in range(4):
#
#     P = obs_names[i]
#     vel_r = doppler_velocity_point(obs_points[i])
#
#     print(f'\t{P}\t{vel_r:10.3f} m/s')
'''
Point	Doppler velocity

	A	 -3744.566 m/s
	B	  1606.148 m/s
	C	  1439.322 m/s
	D	 -2752.272 m/s
'''

'''Creating doppler map'''
# doppler_map('sub FoV', field_of_view='sub', save=True)
# doppler_map('full FoV', save=True)




plt.show()
