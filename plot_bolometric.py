import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as itg
from matplotlib import rcParams
rcParams['font.family'] = 'serif'

# Import Data
phase, lum, lum_error, temp, temp_error, radius, rad_error = np.genfromtxt('Bolometric.txt', unpack = True)
data = phase > 0

# Integrate Luminosity
total_lum=np.log10(itg.trapz(lum[data],phase[data]*3600*24))
total_lum_model=np.log10(itg.trapz(lum[:4],phase[:4]*3600*24))

# Plot Bolometric Lightcurve
plt.figure(figsize=(4,6))
plt.subplots_adjust(hspace=0)
plt.subplot(311)
plt.xlim(-40, 250)
plt.ylim(42.5, 44.6)
plt.errorbar(phase[data], np.log10(lum[data]), yerr = 1/np.log(10) * lum_error[data] / lum[data], fmt = 'o', color = 'g', alpha = 0.4)
plt.fill_between(phase[data], 0, np.log10(lum[data]), color = 'gray', alpha = 0.2, label = r'log($E$ / erg) = %s'%np.around(total_lum, 1), linewidth = 0)
plt.errorbar(phase[:3], np.log10(lum[:3]), yerr = 1/np.log(10) * lum_error[:3] / lum[:3], fmt = '*', color = 'g', alpha = 0.4, ms = 10)
plt.fill_between(phase[:4], 0, np.log10(lum[:4]), color = 'magenta', alpha = 0.2, label = r'log($E$ / erg) = %s'%np.around(total_lum_model, 1), linewidth = 0)
plt.legend(loc = 'upper right', frameon = False, fontsize = 8)
plt.xlabel('Rest days from bolometric maximum')
plt.ylabel(r'log$_{10} (\mathit{L}_{bol}\,/$erg$\,s^{-1})$')
plt.tick_params(axis='both', bottom=False, labelbottom=False)

# Plot Temperature
plt.subplot(312)
plt.xlim(-40, 250)
plt.ylim(12, 32)
plt.errorbar(phase[data], temp[data] / 1000, temp_error[data] / 1000, fmt = 'o', color = 'g', alpha = 0.4)
plt.errorbar(phase[~data], temp[~data] / 1000, temp_error[~data] / 1000, fmt = '*', color = 'g', alpha = 0.4, ms = 10)
plt.ylabel('Temperature / 1000K', labelpad = 12)
plt.tick_params(axis='both', bottom=False, labelbottom=False)

# Plot Radius
plt.subplot(313)
plt.xlim(-40, 250)
plt.ylim(0.0, 1.5)
plt.errorbar(phase[data],  radius[data] / 1E15, rad_error[data] / 1E15, fmt = 'o', color = 'g', alpha = 0.4)
plt.errorbar(phase[~data], radius[~data] / 1E15, rad_error[~data] / 1E15, fmt = '*', color = 'g', alpha = 0.4, label = 'MOSFiT Estimates', ms = 10)
plt.legend(loc = 'upper right', frameon = False, fontsize = 8)
plt.ylabel(r'Radius / $10^{15}$cm')
plt.xlabel('Phase [rest days]')
plt.savefig("Bolometric.pdf", bbox_inches = 'tight')
plt.clf(); plt.close('all')

