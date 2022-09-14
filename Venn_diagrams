import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Import joined data
df = pd.read_csv('data/Ub_SUMO_sites.csv', sep=",",
                      dtype={'PROTEIN': str, 'MOD_RSD': np.int32, 'MW_kD': np.float64, 'ORGANISM': str, '_merge': str})

# Count sites
ub = sum((df['_merge'] == 'Ub'))
sumo = sum((df['_merge'] == 'SUMO'))
both = sum((df['_merge'] == 'both'))

# Group proteins together
df_prot = df.groupby('PROTEIN').agg({'_merge': lambda x: ''.join(x.unique())})

# Count proteins
ub_prot = sum((df_prot['_merge'] == 'Ub'))
sumo_prot = sum((df_prot['_merge'] == 'SUMO'))

# Calculate proteins with SUMO and Ub sites
both_prot = len(df_prot) - (ub_prot + sumo_prot)


# Plot venn diagram of sites
venn2(subsets=(ub, sumo, both), set_labels=['Ub', 'SUMO'])
plt.savefig('Venn_Ub+SUMO_sites.png', dpi=300)

# Plot venn diagram of proteins
#venn2(subsets=(ub_prot, sumo_prot, both_prot), set_labels=['Ub', 'SUMO'])
#plt.savefig('Venn_Ub+SUMO_Protein.png', dpi=300)
