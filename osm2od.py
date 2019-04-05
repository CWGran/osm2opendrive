import numpy as np
import matplotlib.pyplot as plt

from xml.etree.ElementTree import ElementTree
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement

from xml.sax.saxutils import unescape

from geopy import distance
from scipy.optimize import curve_fit

# For debugging
from pprint import pprint

import xml.etree.ElementTree as ET
import argparse

from road import Node, Road


def readNodes(e):
    nodes = {}
    
    # Read all nodes and their coordinates into an array
    for node in e.findall('node'):
        n = Node(node.get("id"), node.get("lat"), node.get("lon"))

        nodes[node.get("id")] = n

    return nodes


def readRoads(e, nodes):
    roads = []

    # Read all roads/ways into an array
    for road in e.findall('way'):
        r = Road(road.get("id"))

        # Read all information about each road
        for tag in road.findall('tag'):
            setattr(r, tag.get("k"), tag.get("v"))
        
        # Connect to the nodes
        for nd in road.findall('nd'):
            r.nodes.append(nodes[nd.get("ref")])

        roads.append(r)

    return roads

def readOSM(filename):
    e = ET.parse('osmmaps/{}'.format(filename)).getroot()

    print("Reading file osmmaps/{}".format(filename))

    nodes = readNodes(e)
    roads = readRoads(e, nodes)

    print("Finished reading file, found {} nodes and {} roads.".format(len(nodes), len(roads)))

    return nodes, roads

def buildXML(filename, roads):

    filename = filename if filename.split(".")[-1] == ".xml" else "".join(filename.split(".")[0]) + ".xml"
    
    print("Building XML output...")

    root = Element('OpenDRIVE')
    tree = ElementTree(root)

    # Setup header record
    header = SubElement(root, "header")
    header.set("revMajor", "1")
    header.set("revMinor", "0")
    header.set("vendor", "Baidu")

    # Setup Geo Reference
    georef = SubElement(header, "geoReference")
    # TODO: Get CDATA working with ElementTree, or switch to lxml.etree
    #georef.text = "<![CDATA[+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs]]>"

    print(len(roads))
    for r in roads:
        pprint(vars(r))

        road = SubElement(root, "road")

        if hasattr(r, "name"):
            road.set("name", r.name)
        else:
            road.set("name", "")

        road.set("id", r.id)

        road.set("junction", "-1")

        rv = SubElement(road, "routeView")
        roadgeo = SubElement(rv, "geometry")
        roadgeo.set("sOffset", "0")
        roadgeo.set("x", str(r.nodes[0].lng))
        roadgeo.set("y", str(r.nodes[0].lat))
        roadgeo.set("z", "0")

        length = 0
        for i in range(len(r.nodes)-1):
            p1 = r.nodes[i]
            p2 = r.nodes[i+1]
            
            length += distance.distance((p1.lat, p1.lng), (p2.lat, p2.lng)).m

        roadgeo.set("length", str(length))

        ps = SubElement(roadgeo, "pointSet")

        # Use this if the first point should not be included
        # for i in range(1, len(r.nodes)):
        for n in r.nodes:
            p = SubElement(ps, "point")
            p.set("x", str(n.lng))
            p.set("y", str(n.lat))
            p.set("z", "0")
        
    print("Done building, writing to '{}'".format(filename))
    tree.write(open('test.xml', 'w'), encoding='unicode')

def find_borders(roads):
    road = 0
    road = roads[road]
    pprint(vars(road))
    points = []

    #for i in range(0, 3):
        #points.append((road.nodes[i].lng, road.nodes[i].lat))

    print("Antall noder: {}".format(len(road.nodes)))
    for n in road.nodes:
        points.append((n.lng, n.lat))

    points = np.array(points)
    vectors = []

    road_width = 5/2
    border_left = []
    border_right = []
    print(points)
    for i in range(len(points)-1):
        p1 = points[i]
        p2 = points[i+1]

        # Vector between the current and the next point
        v = np.array([p2[0]-p1[0], p2[1]-p1[1]])
        # The distance between the two points, calculated in actual meters
        vlen = distance.distance(p1, p2).m

        # The relationship between the meter distance, and the vector norm
        a = vlen/np.linalg.norm(v)
        
        # The two orthogonal vectors, scaled using a and the desired length as specified in road_width
        # Two new points are made by adding the vectors to p1
        left = np.array([-v[1], v[0]])
        l = (road_width/a)*left/np.linalg.norm(left)
        lp = (p1[0] + l[0], p1[1] + l[1])

        right = np.array([v[1], -v[0]])
        r = (road_width/a)*right/np.linalg.norm(right)
        rp = (p1[0] + r[0], p1[1] + r[1])

        border_left.append(lp)
        border_right.append(rp)

        # If this is the last iteration, add a point for the final point using the same orthogonal vectors as for n-1
        if i == len(points)-2:
            lp = (p2[0] + l[0], p2[1] + l[1])
            rp = (p2[0] + r[0], p2[1] + r[1])

            border_left.append(lp)
            border_right.append(rp)

    # A simple print of the calculated meter distance between the old and the generated points
    print("Checking accuracy of generated points:")
    for i in range(len(points)):
        print("Point {}".format(i))
        print("Left: {}".format(distance.distance(points[i], border_left[i]).m))
        print("Right: {}".format(distance.distance(points[i], border_right[i]).m))
        print("Left-right: {}".format(distance.distance(border_left[i], border_right[i]).m))

    plt.plot(points, color="red")
    plt.plot(border_left, color="green")
    plt.plot(border_right, color="blue")
    plt.show()
        

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('file', help="Input filename")

    args = parser.parse_args()

    if args.file:
        filename = args.file

    nodes, roads = readOSM(filename)
    find_borders(roads)
    #buildXML("test", roads)


if __name__ == "__main__":
    main()
