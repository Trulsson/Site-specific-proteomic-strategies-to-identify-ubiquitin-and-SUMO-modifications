import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('Distance_Ub_SUMO.csv', sep=",",
                 usecols=['PROTEIN', 'Ub', 'SUMO'],
                 dtype={'PROTEIN': str, 'Ub': np.int32, 'MW_kD': np.float64, '_merge': str})


# Plot the Figure
labels = ["Ub", "SUMO"]
xlabels = [1]
for i in range(100, 500, 100):
    xlabels.append(i)

plt.rcParams["font.size"] = "26"
plt.rcParams["font.family"] = "Arial"

plt.figure(figsize=(20, 10))
plt.title("Average sequence distance between lysine modifications")

plt.hist(df3, label=labels, bins=50, range=[0, 500])

# X axis limit
ax.set_xlim(1, 30)
ax2.set_xlim(1, 30)

plt.ylabel("Number of proteins")
plt.xlabel("Average distance (AAs)")
plt.legend()
plt.grid(which="major")
plt.savefig("distance_500.png", dpi=300)

#plt.show()
