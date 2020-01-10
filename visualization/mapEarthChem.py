# Place earthchem data points on map, color coded by some feature. 

import h5py
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# from https://stackoverflow.com/questions/4601373/better-way-to-shuffle-two-numpy-arrays-in-unison
def unison_shuffled_copies(a, b, c, d):
    assert len(a) == len(b)
    assert len(b) == len(c)
    p = np.random.permutation(len(a))
    return a[p], b[p], c[p], d[p]

with h5py.File('../igncn1.mat','r') as f:
    lats = f['Latitude'][0]
    longs = f['Longitude'][0]
    #comp = f['Upper_Rho'][0]
    #comp = f['SiO2'][0]
    comp = f['tc1Crust'][0]
    test = f['tc1Crust'][0]


    print(len(lats))

    # Randomly shuffle so overlapping points don't have a trend that biases appearance
    lats, longs, comp, test = unison_shuffled_copies(lats, longs, comp,test)

    #comp = [c*1000 for c in comp]

    lats = [l for i,l in enumerate(lats) if not np.isnan(test[i])]
    longs = [l for i,l in enumerate(longs) if not np.isnan(test[i])]
    comp = [l for i,l in enumerate(comp) if not np.isnan(test[i])]

    print(len(lats))
    # Projection: Cylindrical equal area
    # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
    # are the lat/lon values of the lower left and upper right corners
    # of the map.
    # resolution = 'c' means use crude resolution coastlines.
    plt.figure(figsize=(10, 6), dpi=500)
    m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.drawcoastlines()
    # Draw points
    scatter = m.scatter(longs, lats, latlon=True, marker='o', c=comp, s=2)
    #m.fillcontinents(color='coral',lake_color='aqua')
    m.drawmapboundary(fill_color='grey')

    #plt.legend(*scatter.legend_elements(num=5),title="SiO2 Composition")
    plt.colorbar(orientation="horizontal", label="Depth to 550C Isotherm (km)")

    plt.savefig("../output/tc1-map.png")
