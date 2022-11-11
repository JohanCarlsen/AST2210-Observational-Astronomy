import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

filename = 'ADP.2017-03-27T12_08_50.541.fits'
hdu = fits.open(filename)
# hdu.info()
'''
Filename: ADP.2017-03-27T12_08_50.541.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU    1345   ()
  1  DATA          1 ImageHDU        44   (320, 317, 3682)   float32
  2  STAT          1 ImageHDU        44   (320, 317, 3682)   float32
'''

data = hdu[1].data
hdr = hdu[1].header
# print(data.shape)
'''
(3682, 317, 320)
'''

fluxMean = np.nanmean(data, 0)
# print(fluxMean.shape)
'''
(317, 320)
'''

avg1 = np.nanmean(data, -1)
spect = np.mean(avg1, -1)
print(spect.shape)

plt.figure(figsize=(7, 7))
plt.title(r'Average flux density in the range $\lambda \in [4750,9351] \AA{}$')
im = plt.imshow(np.flip(np.nanmean(data, 0), 0), cmap='gray', vmin=0, vmax=2137)
plt.colorbar(im, fraction=.046, pad=.04, label='Flux density [$10^{-20}$erg s$^{-1}$cm$^{-2}\AA{}^{-1}$]')
plt.tight_layout()
# plt.show()

lambda0 = hdr['CRVAL3']
dLambda = hdr['CD3_3']
lenWave = hdr['NAXIS3']
wavelengths = np.linspace(lambda0, lambda0 + (lenWave - 1) * dLambda, lenWave)

plt.figure(figsize=(8, 4.5))
plt.plot(wavelengths, spect, lw=1, color='black')
plt.show()

lowBound = 6590
upBound = 6610

lowIdx = np.array(np.where(wavelengths >= lowBound))[0, 0]
upIdx = np.array(np.where(wavelengths <= upBound))[-1, -1]

extractedData = data[lowIdx:upIdx]

plt.figure(figsize=(7, 7))
plt.title(f'Average flux density in the range $\lambda \in [{lowBound},{upBound}] \AA$')
im = plt.imshow(np.flip(np.nanmean(extractedData, 0), 0), cmap='gray', vmin=0, vmax=4500)
plt.colorbar(im, fraction=.046, pad=.04, label='Flux density [$10^{-20}$erg s$^{-1}$cm$^{-2}\AA{}^{-1}$]')
plt.tight_layout()
# plt.show()
#
r = 5
centerX = 170
centerY = 150
idx = [centerY - r, centerY + r, centerX -r, centerX + r]

collapsed = np.nansum(data, 0)
apertureData = collapsed[idx[0]:idx[1], idx[2]:idx[3]]
fluxMean_2 = np.mean(apertureData)
