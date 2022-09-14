import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Import joined data
df = pd.read_csv('data/Ub_SUMO_sites.csv', sep=",",
                 usecols=['PROTEIN', 'MOD_RSD', 'MW_kD', '_merge'],
                 dtype={'PROTEIN': str, 'MOD_RSD': np.int32, 'MW_kD': np.float64, '_merge': str})

# Group proteins and count number of sites
df_Ub = df[df['_merge'] != 'SUMO'].groupby('PROTEIN')['_merge'].count()
df_SUMO = df[df['_merge'] != 'Ub'].groupby('PROTEIN')['_merge'].count()

labels = ["Ub", "SUMO"]

plt.rcParams['font.size'] = '30'
plt.rcParams['font.family'] = 'Arial'


# Create 2 subplots and set size of upper plot (outliers) and lower plot (most of the data)
fig, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 5]})

bins = int((max(df_Ub)))

# Plot the data
ax.hist([df_Ub, df_SUMO], bins=bins, label=labels, rwidth=0.8)
ax2.hist([df_Ub, df_SUMO], bins=bins, label=labels, rwidth=0.8)


# Divide what data to show on each subplot
ax.set_ylim(2000, 3000)
ax2.set_ylim(0, 1500)

# X axis limit
ax.set_xlim(1, 30)
ax2.set_xlim(1, 30)

# Hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # Don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .010  # How big to make the diagonal lines in axes coordinates
# Arguments to pass to plot diagonal lines on y-axes
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
ax.grid(which='major')

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
ax2.grid(which='major')

fig.set_figheight(10)
fig.set_figwidth(30)

plt.xlabel("Modified lysines")
plt.ylabel("Number of proteins")
plt.legend()
#plt.show()
plt.savefig('phosphosite_histogram_brokenaxis.png', dpi=300)
