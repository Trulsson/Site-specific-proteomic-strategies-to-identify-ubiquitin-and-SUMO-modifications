import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

df_SUMO = pd.read_csv('data/Sumoylation_site_dataset.txt', sep="\t",
                      usecols=['PROTEIN', 'MOD_RSD', 'MW_kD', 'ORGANISM'],
                      dtype={'PROTEIN': str, 'MOD_RSD': str, 'MW_kD': np.float64, 'ORGANISM': str})
df_Ub = pd.read_csv('data/Ubiquitination_site_dataset.txt', sep="\t",
                    usecols=['PROTEIN', 'MOD_RSD', 'MW_kD', 'ORGANISM'],
                    dtype={'PROTEIN': str, 'MOD_RSD': str, 'MW_kD': np.float64, 'ORGANISM': str})


df_SUMO["MOD_RSD"].replace('K', '', regex=True, inplace=True)
df_SUMO.replace('(-[s,m,u,b]{2})', '', regex=True, inplace=True)
df_SUMO_human = df_SUMO.loc[df_SUMO["ORGANISM"] == "human"]

df_Ub_k = df_Ub[df_Ub['MOD_RSD'].str.startswith("K")]
df_Ub_k["MOD_RSD"].replace('K', '', regex=True, inplace=True)
df_Ub_k.replace('(-[s,m,u,b]{2})', '', regex=True, inplace=True)
df_Ub_human = df_Ub_k.loc[df_Ub["ORGANISM"] == "human"]

df_outer = pd.merge(df_Ub_human, df_SUMO_human, how='outer',
                    on=['PROTEIN', 'MOD_RSD', 'ORGANISM', 'MW_kD'], indicator=True)
df_inner = pd.merge(df_Ub_human, df_SUMO_human, how='inner',
                    on=['PROTEIN', 'MOD_RSD', 'ORGANISM', 'MW_kD'], indicator=True)
df_both = pd.concat([df_inner, df_outer], ignore_index=True)
df_both.replace('left_only', 'Ub', inplace=True)
df_both.replace('right_only', 'SUMO', inplace=True)

df_both.to_csv('data/Ub_SUMO_sites.csv')
