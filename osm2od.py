import numpy as np
import matplotlib.pyplot as plt
import math

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

    for r in roads:
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

        roadgeo.set("length", str(road_length(r.nodes)))

        ps = SubElement(roadgeo, "pointSet")

        # Not sure why these points are here as the reference line is supposed to be stored within <center>
        # Need to look at an example
        # Use this if the first point should not be included
        # for i in range(1, len(r.nodes)):
        for n in r.nodes:
            p = SubElement(ps, "point")
            p.set("x", str(n.lng))
            p.set("y", str(n.lat))
            p.set("z", "0")

        ##########
        #### Lanes
        ##########

        lanes = SubElement(road, "lanes")

        num_lanes = 1
        lane_width = 3.0
        if hasattr(r, "lanes"):
            num_lanes = int(r.lanes)

        laneSec = SubElement(lanes, "laneSection")
        laneSec.set("singleSide", "true")          # True if both directions share the same laneSection

        boundaries = SubElement(laneSec, "boundaries")

        # If the lane number is odd and greater than 1, only care about num_lanes-1 lanes
        if num_lanes > 1 and num_lanes % 2 != 0:
            num_lanes -= 1

        # Lane boundaries
        left_boundary = SubElement(boundaries, "boundary")
        left_boundary.set("type", "leftBoundary")

        right_boundary = SubElement(boundaries, "boundary")
        right_boundary.set("type", "rightBoundary")

        # Lane boundary geometries
        leftb_geo = SubElement(left_boundary, "geometry")
        leftb_geo_ps = SubElement(leftb_geo, "pointSet")

        rightb_geo = SubElement(right_boundary, "geometry")
        rightb_geo_ps = SubElement(rightb_geo, "pointSet")
        
        boundary_width = num_lanes/2.0*lane_width
        left_boundary_points = find_parallel(r.nodes, boundary_width, True)
        right_boundary_points = find_parallel(r.nodes, boundary_width, False)

        for i in range(len(r.nodes)):
            # Left
            lp = SubElement(leftb_geo_ps, "point")
            lp.set("x", str(left_boundary_points[i][0]))
            lp.set("y", str(left_boundary_points[i][1]))
            lp.set("z", "0")

            # Right
            rp = SubElement(rightb_geo_ps, "point")
            rp.set("x", str(right_boundary_points[i][0]))
            rp.set("y", str(right_boundary_points[i][1]))
            rp.set("z", "0")

        # Center is supposed to store the reference line
        # Left/right stores the borders of left/right lanes

        center = SubElement(laneSec, "center")
        center_lane = SubElement(center, "lane")

        center_lane.set("id", "0")
        center_lane.set("uid", "{}_0".format(r.id))
        center_lane.set("type", "none")
        center_lane.set("direction", "bidirection")
        center_lane.set("turnType", "noTurn")   # Not sure what this means

        center_line = SubElement(center_lane, "centerLine")
        cl_geo = SubElement(center_line, "geometry")
        cl_geo.set("sOffset", "0")
        cl_geo.set("x", str(r.nodes[0].lng))
        cl_geo.set("y", str(r.nodes[0].lat))
        cl_geo.set("z", "0")
        cl_geo.set("length", str(road_length(r.nodes)))

        cl_geo_ps = SubElement(cl_geo, "pointSet")

        for n in r.nodes:
            p = SubElement(cl_geo_ps, "point")
            p.set("x", str(n.lng))
            p.set("y", str(n.lat))
            p.set("z", "0")

        right = SubElement(laneSec, "right")

        if num_lanes > 1:
            left = SubElement(laneSec, "left")

        for i in range(math.ceil(num_lanes/2)):
            # Right, only add this if num_lanes == 1
            right_lane = SubElement(right, "lane")
            right_lane.set("id", "{}".format(i+1))
            right_lane.set("uid", "{}_0{}".format(r.id, i+1))
            right_lane.set("type", "driving")
            right_lane.set("direction", "bidirection" if num_lanes == 1 else "forward")
            right_lane.set("turnType", "noTurn")    # Not sure what this means

            right_border = SubElement(right_lane, "border")
            right_border.set("virtual", "TRUE")     # "Identify whether the lane boundary exists in real world"

            right_border_points = find_parallel(r.nodes, i*lane_width, True)

            rb_geo = SubElement(right_border, "geometry")
            rb_geo.set("sOffset", "0")
            rb_geo.set("x", str(right_border_points[0][0]))
            rb_geo.set("y", str(right_border_points[0][1]))
            rb_geo.set("z", "0")
            rb_geo.set("length", str(road_length(right_border_points)))

            rb_geo_ps = SubElement(rb_geo, "pointSet")

            for n in right_border_points:
                p = SubElement(rb_geo_ps, "point")
                p.set("x", str(n[0]))
                p.set("y", str(n[1]))
                p.set("z", "0")

            if num_lanes > 1:
                left_lane = SubElement(left, "lane")
                left_lane.set("id", "-{}".format(i+1))
                left_lane.set("uid", "{}_1{}".format(r.id, i+1))
                left_lane.set("type", "driving")
                left_lane.set("direction", "backward")
                left_lane.set("turnType", "noTurn")    # Not sure what this means

                left_border = SubElement(left_lane, "border")
                left_border.set("virtual", "TRUE")     # "Identify whether the lane boundary exists in real world"

                left_border_points = find_parallel(r.nodes, i*lane_width, False)

                lb_geo = SubElement(left_border, "geometry")
                lb_geo.set("sOffset", "0")
                lb_geo.set("x", str(left_border_points[0][0]))
                lb_geo.set("y", str(left_border_points[0][1]))
                lb_geo.set("z", "0")
                lb_geo.set("length", str(road_length(left_border_points)))

                lb_geo_ps = SubElement(lb_geo, "pointSet")

                for n in left_border_points:
                    p = SubElement(lb_geo_ps, "point")
                    p.set("x", str(n[0]))
                    p.set("y", str(n[1]))
                    p.set("z", "0")
        
    print("Done building, writing to '{}'".format(filename))
    tree.write(open('test.xml', 'w'), encoding='unicode')

# Calculate road length
def road_length(road):
    length = 0
    for i in range(len(road)-1):
        p1 = road[i]
        p2 = road[i+1]
        
        if isinstance(p1, Node):
            length += distance.distance((p1.lat, p1.lng), (p2.lat, p2.lng)).m
        else:
            length += distance.distance(p1, p2).m

    return length
    
def find_parallel(road, width, left):
    points = []

    #for i in range(0, 3):
        #points.append((road.nodes[i].lng, road.nodes[i].lat))

    for n in road:
        points.append((n.lng, n.lat))

    points = np.array(points)
    vectors = []

    road_width = 5/2
    parallel = []
    for i in range(len(points)-1):
        p1 = points[i]
        p2 = points[i+1]

        # Vector between the current and the next point
        v = np.array([p2[0]-p1[0], p2[1]-p1[1]])
        # The distance between the two points, calculated in actual meters
        vlen = distance.distance(p1, p2).m

        # The relationship between the meter distance, and the vector norm
        a = vlen/np.linalg.norm(v)
        
        # The two orthogonal vectors, scaled using a and the desired length as specified in width
        # Two new points are made by adding the vectors to p1
        lv = np.array([-v[1], v[0]]) if left else np.array([v[1], -v[0]])
        l = (width/a)*lv/np.linalg.norm(lv)
        lp = (p1[0] + l[0], p1[1] + l[1])
        
        parallel.append(lp)

        # If this is the last iteration, add a point for the final point using the same orthogonal vectors as for n-1
        if i == len(points)-2:
            lp = (p2[0] + l[0], p2[1] + l[1])

            parallel.append(lp)
    
    return parallel

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('file', help="Input filename")

    args = parser.parse_args()

    if args.file:
        filename = args.file

    nodes, roads = readOSM(filename)
    #find_parallel(roads)
    buildXML("test", roads)


if __name__ == "__main__":
    main()
