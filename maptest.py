import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree
from scipy.interpolate import interp1d, CubicSpline
from mpl_toolkits.basemap import Basemap

nodes = {}
lats = []
lons = []

e = xml.etree.ElementTree.parse('map2').getroot()

for node in e.findall('node'):
    nodes[node.get("id")] = [node.get("lat"), node.get("lon")]
    lats.append(float(node.get("lat")))
    lons.append(float(node.get("lon")))

ways = []
for way in e.findall('way'):
    lat = []
    lon = []

    for nd in way.findall('nd'):
        node = nodes[nd.get("ref")]
        lat.append(node[0])
        lon.append(node[1])

    ways.append([lat, lon])

print("Found {} ways and {} nodes".format(len(ways), len(nodes)))

my_map = Basemap(projection='merc', lat_0=63.45, lon_0=10.42, resolution='h', area_thresh=0.1, llcrnrlon=10.43, llcrnrlat=63.4, urcrnrlon=10.2, urcrnrlat=63.48)
my_map.bluemarble()
my_map.drawcoastlines()
my_map.drawcountries()
my_map.drawmapboundary()

x, y = my_map(lons, lats)
my_map.plot(x, y, 'ro', markersize=1)

plt.show()
