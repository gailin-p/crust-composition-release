import h5py
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

with h5py.File('../igncn1.mat','r') as f:
    lats = f['Latitude'][0]
    longs = f['Longitude'][0]

    # Projection: Cylindrical equal area
    # llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
    # are the lat/lon values of the lower left and upper right corners
    # of the map.
    # resolution = 'c' means use crude resolution coastlines.
    plt.figure(figsize=(20, 12))
    m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,\
                llcrnrlon=-180,urcrnrlon=180,resolution='c')
    m.drawcoastlines()
    # Draw points
    m.scatter(longs, lats, latlon=True, marker='o', alpha=.2)
    #m.fillcontinents(color='coral',lake_color='aqua')
    m.drawmapboundary(fill_color='aqua')
    plt.title("EarthChem data locations")
    plt.show()
