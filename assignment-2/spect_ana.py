# This program will take the spectrum file "IRAS13120_spec.txt" and fit
# a Gaussian to it.
import numpy as np
import matplotlib.pyplot as plt

_v, antennaT = np.loadtxt('IRAS13120_spec.txt', unpack=True)
v = _v / 1e3

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(v, antennaT, color='black' , lw=.5)

ax.set_xlabel('Rad. vel. [$m/s$]')
plt.show()
