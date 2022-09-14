import pandas as pd
import numpy as np


def check_distance(a_list, x):
    # Takes a list of numbers (a_list) and a number (x) and returns the distance of the number to each value in the list

    distance = []
    for i, number in enumerate(a_list):
        if int(number) > int(x):
            dist = int(number) - int(x)
            # print(f"The number in the list {number} - the input value {x} = {dist}")
        else:
            dist = int(x) - int(number)
            # print(f"The input value {x} - the number in the list {number} = {dist}")
        distance.append(dist)
    result = sum(distance) / len(distance)
    # print(f"The average of {distance} is {result}")
    return result


# Import joined data
df = pd.read_csv('data/Ub_SUMO_sites.csv', sep=",",
                 usecols=['PROTEIN', 'MOD_RSD', 'MW_kD', '_merge'],
                 dtype={'PROTEIN': str, 'MOD_RSD': np.int32, 'MW_kD': np.float64, '_merge': str})

# Separate Ub and SUMO sites
df_Ub = df[df['_merge'] != 'SUMO']
df_SUMO = df[df['_merge'] != 'Ub']

# Group all modifications of the same proteins together
Ub_grouped = df_Ub.groupby(by=["PROTEIN"])
SUMO_grouped = df_SUMO.groupby(by=["PROTEIN"])

avg_distance_dict_Ub = {}
avg_distance_dict_SUMO = {}

# Calculate distance between modifications within the same gene for Ub
for key, item in Ub_grouped:
    avg_distance = []
    lysine_pos = []
    # Take only proteins with more than 1 modification
    if len(item) > 1:
        # Iterate over proteins and extract lysine positions
        for row in Ub_grouped["MOD_RSD"].get_group(key):
            lysine_pos.append(row)
            # Calculate the distance of each lysine position from all modified lysine positions of that protein
            for idx, position in enumerate(lysine_pos):
                avg_distance.append(check_distance(lysine_pos, position))
        # Calculate average distance within protein
        avg_distance_protein = sum(avg_distance) / len(avg_distance)
        avg_distance_dict_Ub.update({key: avg_distance_protein})
    else:
        pass

# Calculate distance between modifications within the same gene for SUMO
for key, item in SUMO_grouped:
    avg_distance = []
    lysine_pos = []
    # Take only proteins with more than 1 modification
    if len(item) > 1:
        # Iterate over proteins and extract lysine positions
        for row in SUMO_grouped["MOD_RSD"].get_group(key):
            lysine_pos.append(row)
            # Calculate the distance of each lysine position from all modified lysine positions of that protein
            for idx, position in enumerate(lysine_pos):
                avg_distance.append(check_distance(lysine_pos, position))
        # Calculate average distance within protein
        avg_distance_protein = sum(avg_distance) / len(avg_distance)
        avg_distance_dict_SUMO.update({key: avg_distance_protein})
    else:
        pass


# Put all the average distances in the same dataframe
df_Ub_distance = pd.DataFrame.from_dict(avg_distance_dict_Ub, orient='index', columns=['Ub'])
df_Ub_distance = df_Ub_distance.rename_axis('PROTEIN').reset_index()

df_SUMO_distance = pd.DataFrame.from_dict(avg_distance_dict_SUMO, orient='index', columns=['SUMO'])
df_SUMO_distance = df_SUMO_distance.rename_axis('PROTEIN').reset_index()

df_outer = pd.merge(df_Ub_distance, df_SUMO_distance, how='outer', on='PROTEIN')
df_inner = pd.merge(df_Ub_distance, df_SUMO_distance, how='inner', on='PROTEIN')
df_both = pd.concat([df_inner, df_outer], ignore_index=True)

print(df_both)
df_both.to_csv('Distance_Ub_SUMO.csv')
