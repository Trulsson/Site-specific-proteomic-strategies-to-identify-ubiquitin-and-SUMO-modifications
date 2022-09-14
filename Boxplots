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

print(df_Ub.describe())
print(df_SUMO.describe())

labels = ["Ub", "SUMO"]

# Set font and font size
plt.rcParams['font.size'] = '22'
plt.rcParams['font.family'] = 'Arial'

# Create 2 subplots and set size of upper plot (outliers) and lower plot (most of the data)
f, (ax, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 5]})

# Plot the data
ax.boxplot([df_Ub, df_SUMO], labels=labels)
ax2.boxplot([df_Ub, df_SUMO], labels=labels)

# Divide what data to show on each subplot
ax.set_ylim(480, 500)  # outliers only
ax2.set_ylim(0, 220)  # most of the data
y_ticks = [1]
for i in range(30, 220, 30):
    y_ticks.append(i)
ax2.set(yticks=y_ticks)

# Hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # Don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .015  # How big to make the diagonal lines in axes coordinates
# Arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

plt.savefig("phosphosite_boxplot.png", dpi=300)
