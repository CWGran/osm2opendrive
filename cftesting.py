import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree
from scipy.interpolate import interp1d, CubicSpline, splprep, splev
from pyproj import Proj, transform

def convert_coords(lats, lons):
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init='epsg:3857')
    x, y = transform(inProj, outProj, lats, lons)
    return x, y

osmfile = "osmmaps/" + "midtbyen"

unsupported = ["pedestrian", "unclassified", "cycleway", "footway", "path"]
driveable = ["motorway", "trunk", "primary", "secondary", "tertiary", "residential", "service", "living_street", "track", "road", "unclassified"]
for r in driveable.copy():
    driveable.append(r + "_link")

nodes = {}
lats = []
lons = []

e = xml.etree.ElementTree.parse(osmfile).getroot()

for node in e.findall('node'):
    nodes[node.get("id")] = [node.get("lon"), node.get("lat")]
    lats.append(float(node.get("lat")))
    lons.append(float(node.get("lon")))

ways = []
for way in e.findall('way'):
    lat = []
    lon = []

    supp = False
    name = False
    for tag in way.findall('tag'):
        if tag.get('k') == "highway":
            if tag.get('v') in driveable:
                supp = True
        elif tag.get('k') == "name":
            name = tag.get('v')
            print("Way found: {}".format(name))
    if not supp or not name:
        continue

    for nd in way.findall('nd'):
        node = nodes[nd.get("ref")]
        lat.append(float(node[0]))
        lon.append(float(node[1]))

    x, y = convert_coords(np.array(lat), np.array(lon))

    x0 = x[0]
    y0 = y[0]

    # Dette blir for enkelt tror jeg
    #x = x-x0
    #y = y-y0


    #ways.append(np.dstack((x,y)))
    ways.append(np.array(list(zip(x,y))))


for w in ways:
    k = 3
    if len(w) <= 3:
        k = len(w)-1

    tck, u = splprep(w.T, u=None, s=0.0, k=k) 
    u_new = np.linspace(u.min(), u.max(), 1000)
    x_new, y_new = splev(u_new, tck, der=0)

    plt.plot(w[:,0], w[:,1], 'o')
    plt.plot(x_new, y_new, 'b--')

    #plt.scatter(w[:,0], w[:,1])

print("Found {} ways and {} nodes".format(len(ways), len(nodes)))
plt.show()
