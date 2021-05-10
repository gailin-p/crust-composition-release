import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import stats 

nbdat = pd.read_csv("data/remote/latlong_nobin/output/results-upper-meanvprhorange-earthchem-Crust1.0.csv")
bdat = pd.read_csv("data/remote/latlong_weighted/output/results-upper-meanvprhorange-earthchem-Crust1.0.csv")

print(bdat.label)

fig, axs = plt.subplots(2)

axs[0].plot(nbdat["sample_geotherm"], nbdat["SiO2"], 'b.', zorder=1)
nb_means, nb_edges, nbnum = stats.binned_statistic(nbdat["sample_geotherm"], nbdat["SiO2"])
axs[0].hlines(nb_means, nb_edges[:-1], nb_edges[1:], colors='g', lw=5, zorder=2)

axs[1].plot(bdat["sample_geotherm"], bdat["SiO2"], 'b.', zorder=1)
b_means, b_edges, bnum = stats.binned_statistic(bdat["sample_geotherm"], bdat["SiO2"])
axs[1].hlines(b_means, b_edges[:-1], bg_edges[1:], colors='g', lw=5, zorder=2)

