import matplotlib.pyplot as plt
import numpy as np
from astropy import table
from matplotlib.pyplot import cm 
from matplotlib import rcParams
rcParams['font.family'] = 'serif'

# Import Model
model_MJD = np.genfromtxt('models/MJD.txt')

UVM2_Model = np.genfromtxt('models/UVM2_Mag.txt')
B_Model    = np.genfromtxt('models/B_Mag.txt')
g_Model    = np.genfromtxt('models/g_Mag.txt')
i_Model    = np.genfromtxt('models/i_Mag.txt')
r_Model    = np.genfromtxt('models/r_Mag.txt')
U_Model    = np.genfromtxt('models/U_Mag.txt')
UVW1_Model = np.genfromtxt('models/UVW1_Mag.txt')
UVW2_Model = np.genfromtxt('models/UVW2_Mag.txt')
V_Model    = np.genfromtxt('models/V_Mag.txt')

# Phase
MJD0 = 58429
redshift = 0.045726
Phase = (model_MJD - MJD0) / (1 + redshift)

# Import Data
data_table = table.Table.read('Photometry.txt', format = 'ascii')

# Band Offsets and Colors
offsets = np.arange(0, 1.3 * 9, 1.3) 
colors = cm.rainbow(np.linspace(0,1,9))
color_array = {'i':'maroon','r':'r','V':colors[5],'g':'g','B':colors[4],'U':colors[3],'UVW1':colors[2],'UVM2':colors[1],'UVW2':colors[0]}

# Plot
plt.gca().invert_yaxis()
plt.gca().set_xlim(-50, 600)
#plt.gca().set_ylim(bottom=25, top=17)
plt.gca().set_ylim(bottom=31.9, top=15.5)
plt.gca().set_xticks(np.arange(-50, 600, 50))
plt.gca().set_yticks(np.arange(15.5, 31.9, 1))
plt.gca().set_xlabel('Phase [rest days]')
plt.gca().set_ylabel('Apparent Magnitude + Constant')

plt.plot(Phase, i_Model    + offsets[0], color = 'maroon' , alpha=0.3, linewidth=0.2)
plt.plot(Phase, r_Model    + offsets[1], color = 'r'      , alpha=0.3, linewidth=0.2)
plt.plot(Phase, V_Model    + offsets[2], color = colors[5], alpha=0.3, linewidth=0.2)
plt.plot(Phase, g_Model    + offsets[3], color = 'g'      , alpha=0.3, linewidth=0.2)
plt.plot(Phase, B_Model    + offsets[4], color = colors[4], alpha=0.3, linewidth=0.2)
plt.plot(Phase, U_Model    + offsets[5], color = colors[3], alpha=0.3, linewidth=0.2)
plt.plot(Phase, UVW1_Model + offsets[6], color = colors[2], alpha=0.3, linewidth=0.2)
plt.plot(Phase, UVM2_Model + offsets[7], color = colors[1], alpha=0.3, linewidth=0.2)
plt.plot(Phase, UVW2_Model + offsets[8], color = colors[0], alpha=0.3, linewidth=0.2)

# Plot Data
for i, band in enumerate(color_array):
    match = data_table['Filter'] == band
    plt.errorbar((data_table['MJD'][match]- MJD0) / (1 + redshift), data_table['AB_Mag'][match]+offsets[i], yerr=data_table['Error'][match], color=color_array[band], fmt='o', label = band,
                 markeredgecolor='black', markeredgewidth=1, capsize=1,elinewidth=1.5, capthick=2, zorder=10, ms = 10)
    plt.errorbar((data_table['MJD'][match]- MJD0) / (1 + redshift), data_table['AB_Mag'][match]+offsets[i], yerr=data_table['Error'][match], color='k', fmt='o', capsize=2,
                 elinewidth=2.5, capthick=3, zorder=5, ms = 10)

plt.legend(ncol=4, bbox_to_anchor=(-0.12, 1.02),loc='lower left', frameon = True)
plt.savefig('lc_mosfit.pdf', bbox_inches = 'tight')
plt.clf()

