import numpy as np
import matplotlib.pyplot as plt
import math
import datetime
import argparse

from geopy import distance
from scipy.optimize import curve_fit
from lxml import etree

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
    driveable = ["motorway", "trunk", "primary", "secondary", "tertiary", "residential", "service", "living_street", "track", "road", "unclassified"]
    for r in driveable.copy():
        driveable.append(r + "_link")

    # Read all roads/ways into an array
    for road in e.findall('way'):
        r = Road(road.get("id"))

        supp = False

        # Read all information about each road
        for tag in road.findall('tag'):
            setattr(r, tag.get("k"), tag.get("v"))

            if tag.get('k') == "highway":
                if tag.get('v') in driveable:
                    supp = True
        if not supp:
            continue
        
        # Connect to the nodes
        for nd in road.findall('nd'):
            r.nodes.append(nodes[nd.get("ref")])

        roads.append(r)

    return roads

def readOSM(filename):
    e = etree.parse('osmmaps/{}'.format(filename)).getroot()

    print("Reading file osmmaps/{}".format(filename))

    nodes = readNodes(e)
    roads = readRoads(e, nodes)

    print("Finished reading file, found {} nodes and {} roads.".format(len(nodes), len(roads)))

    return nodes, roads

def format_coord(n):
    return "{:.9e}".format(n)

def buildXML(filename, roads, pretty):

    name = filename.split(".")[0]
    #filename = name + ".xml"
    filename = "base_map.xml"
    
    print("Building XML output...")

    root = etree.Element('OpenDRIVE')
    root.set("xmlns", "http://www.opendrive.org")
    tree = etree.ElementTree(root)

    # Setup header record
    header = etree.SubElement(root, "header")
    header.set("revMajor", "1")
    header.set("revMinor", "0")
    header.set("vendor", "Baidu")
    header.set("name", name)
    header.set("version", "1.0")
    header.set("date", datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))

    # Maximum and minimum coordinate values
    # North, south, east, west
    max_coord = [None,None,None,None]

    # Setup Geo Reference
    georef = etree.SubElement(header, "geoReference")
    # TODO: Get CDATA working with ElementTree, or switch to lxml.etree
    georef.text = etree.CDATA("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

    for r in roads:
        road = etree.SubElement(root, "road")

        if hasattr(r, "name"):
            road.set("name", r.name)
        else:
            road.set("name", "")

        road.set("id", r.id)

        road.set("junction", "-1")

        # Routeview
        # This is optional
        #rv = SubElement(road, "routeView")
        #roadgeo = SubElement(rv, "geometry")
        #roadgeo.set("sOffset", "0")
        #roadgeo.set("x", str(r.nodes[0].lng))
        #roadgeo.set("y", str(r.nodes[0].lat))
        #roadgeo.set("z", "0")

        #roadgeo.set("length", str(road_length(r.nodes)))

        #ps = SubElement(roadgeo, "pointSet")

        ## Not sure why these points are here as the reference line is supposed to be stored within <center>
        ## Need to look at an example
        ## Use this if the first point should not be included
        ## for i in range(1, len(r.nodes)):
        #for n in r.nodes:
            #p = SubElement(ps, "point")
            #p.set("x", str(n.lng))
            #p.set("y", str(n.lat))
            #p.set("z", "0")

        ##########
        #### Lanes
        ##########

        lanes = etree.SubElement(road, "lanes")

        num_lanes = 1
        lane_width = 3.0
        if hasattr(r, "lanes"):
            num_lanes = int(r.lanes)

        laneSec = etree.SubElement(lanes, "laneSection")
        laneSec.set("singleSide", "true")          # True if both directions share the same laneSection

        boundaries = etree.SubElement(laneSec, "boundaries")

        # If the lane number is odd and greater than 1, only care about num_lanes-1 lanes
        if num_lanes > 1 and num_lanes % 2 != 0:
            num_lanes -= 1

        # Lane boundaries
        left_boundary = etree.SubElement(boundaries, "boundary")
        left_boundary.set("type", "leftBoundary")

        right_boundary = etree.SubElement(boundaries, "boundary")
        right_boundary.set("type", "rightBoundary")

        # Lane boundary geometries
        leftb_geo = etree.SubElement(left_boundary, "geometry")
        leftb_geo_ps = etree.SubElement(leftb_geo, "pointSet")

        rightb_geo = etree.SubElement(right_boundary, "geometry")
        rightb_geo_ps = etree.SubElement(rightb_geo, "pointSet")
        
        boundary_width = num_lanes/2.0*lane_width
        left_boundary_points = find_parallel(r.nodes, boundary_width, True)
        right_boundary_points = find_parallel(r.nodes, boundary_width, False)

        for i in range(len(r.nodes)):
            # Left
            lp = etree.SubElement(leftb_geo_ps, "point")
            lp.set("x", format_coord(left_boundary_points[i][0]))
            lp.set("y", format_coord(left_boundary_points[i][1]))
            lp.set("z", format_coord(0.0))

            # Right
            rp = etree.SubElement(rightb_geo_ps, "point")
            rp.set("x", format_coord(right_boundary_points[i][0]))
            rp.set("y", format_coord(right_boundary_points[i][1]))
            rp.set("z", format_coord(0.0))

        # Center is supposed to store the reference line
        # Left/right stores the borders of left/right lanes

        center = etree.SubElement(laneSec, "center")
        center_lane = etree.SubElement(center, "lane")

        center_lane.set("id", "0")
        center_lane.set("uid", "{}_0".format(r.id))
        center_lane.set("type", "none")
        #center_lane.set("direction", "bidirection")
        #center_lane.set("turnType", "noTurn")   # Not sure what this means

        center_border = etree.SubElement(center_lane, "border")
        center_border.set("virtual", "TRUE")
        cl_geo = etree.SubElement(center_border, "geometry")
        cl_geo.set("sOffset", "0")
        cl_geo.set("x", format_coord(r.nodes[0].lng))
        cl_geo.set("y", format_coord(r.nodes[0].lat))
        cl_geo.set("z", format_coord(0.0))
        cl_geo.set("length", str(road_length(r.nodes)))

        cl_geo_ps = etree.SubElement(cl_geo, "pointSet")

        for n in r.nodes:
            p = etree.SubElement(cl_geo_ps, "point")
            p.set("x", format_coord(n.lng))
            p.set("y", format_coord(n.lat))
            p.set("z", format_coord(0.0))

            # Check for min/max values:
            # North
            if max_coord[0] == None or max_coord[0] < n.lat:
                max_coord[0] = n.lat

            # South
            if max_coord[1] == None or max_coord[1] > n.lat:
                max_coord[1] = n.lat

            # East
            if max_coord[2] == None or max_coord[2] < n.lng:
                max_coord[2] = n.lng

            # West
            if max_coord[3] == None or max_coord[3] > n.lng:
                max_coord[3] = n.lng

        right = etree.SubElement(laneSec, "right")

        if num_lanes > 1:
            left = etree.SubElement(laneSec, "left")

        for i in range(math.ceil(num_lanes/2)):
            # Right, only add this if num_lanes == 1
            right_lane = etree.SubElement(right, "lane")
            right_lane.set("id", "-{}".format(i+1))
            right_lane.set("uid", "{}_1{}".format(r.id, i+1))
            right_lane.set("type", "driving")
            right_lane.set("direction", "bidirection" if num_lanes == 1 else "forward")
            right_lane.set("turnType", "noTurn")    # Not sure what this means

            # Lane center
            right_center = etree.SubElement(right_lane, "centerLine")
                
            center_pos = (i-1)*lane_width+(lane_width/2)
            right_center_points = find_parallel(r.nodes, center_pos, True)

            rc_geo = etree.SubElement(right_center, "geometry")
            rc_geo.set("sOffset", "0")
            rc_geo.set("x", format_coord(right_center_points[0][0]))
            rc_geo.set("y", format_coord(right_center_points[0][1]))
            rc_geo.set("z", format_coord(0.0))
            rc_geo.set("length", str(road_length(right_center_points)))

            rc_geo_ps = etree.SubElement(rc_geo, "pointSet")

            for n in right_center_points:
                p = etree.SubElement(rc_geo_ps, "point")
                p.set("x", format_coord(n[0]))
                p.set("y", format_coord(n[1]))
                p.set("z", format_coord(0.0))

            # Lane border
            right_border = etree.SubElement(right_lane, "border")
            right_border.set("virtual", "TRUE")     # "Identify whether the lane boundary exists in real world"

            right_border_points = find_parallel(r.nodes, i*lane_width, True)

            rb_geo = etree.SubElement(right_border, "geometry")
            rb_geo.set("sOffset", "0")
            rb_geo.set("x", format_coord(right_border_points[0][0]))
            rb_geo.set("y", format_coord(right_border_points[0][1]))
            rb_geo.set("z", format_coord(0.0))
            rb_geo.set("length", str(road_length(right_border_points)))

            rb_geo_ps = etree.SubElement(rb_geo, "pointSet")

            for n in right_border_points:
                p = etree.SubElement(rb_geo_ps, "point")
                p.set("x", format_coord(n[0]))
                p.set("y", format_coord(n[1]))
                p.set("z", format_coord(0.0))

            if num_lanes > 1:
                left_lane = etree.SubElement(left, "lane")
                left_lane.set("id", "{}".format(i+1))
                left_lane.set("uid", "{}_0{}".format(r.id, i+1))
                left_lane.set("type", "driving")
                left_lane.set("direction", "backward")
                left_lane.set("turnType", "noTurn")    # Not sure what this means

                # Lane center
                left_center = etree.SubElement(left_lane, "centerLine")
                    
                center_pos = (i-1)*lane_width+(lane_width/2)
                left_center_points = find_parallel(r.nodes, center_pos, False)

                lc_geo = etree.SubElement(left_center, "geometry")
                lc_geo.set("sOffset", "0")
                lc_geo.set("x", format_coord(left_center_points[0][0]))
                lc_geo.set("y", format_coord(left_center_points[0][1]))
                lc_geo.set("z", format_coord(0.0))
                lc_geo.set("length", str(road_length(left_center_points)))

                lc_geo_ps = etree.SubElement(lc_geo, "pointSet")

                for n in left_center_points:
                    p = etree.SubElement(lc_geo_ps, "point")
                    p.set("x", format_coord(n[0]))
                    p.set("y", format_coord(n[1]))
                    p.set("z", format_coord(0.0))

                # Lane border
                left_border = etree.SubElement(left_lane, "border")
                left_border.set("virtual", "TRUE")     # "Identify whether the lane boundary exists in real world"

                left_border_points = find_parallel(r.nodes, i*lane_width, False)

                lb_geo = etree.SubElement(left_border, "geometry")
                lb_geo.set("sOffset", "0")
                lb_geo.set("x", format_coord(left_border_points[0][0]))
                lb_geo.set("y", format_coord(left_border_points[0][1]))
                lb_geo.set("z", format_coord(0.0))
                lb_geo.set("length", str(road_length(left_border_points)))

                lb_geo_ps = etree.SubElement(lb_geo, "pointSet")

                for n in left_border_points:
                    p = etree.SubElement(lb_geo_ps, "point")
                    p.set("x", format_coord(n[0]))
                    p.set("y", format_coord(n[1]))
                    p.set("z", format_coord(0.0))
        
    header.set("north", format_coord(max_coord[0]))
    header.set("south", format_coord(max_coord[1]))
    header.set("east", format_coord(max_coord[2]))
    header.set("west", format_coord(max_coord[3]))
    print("Done building, writing to '{}'".format(filename))

    tree.write(filename, xml_declaration=True, pretty_print=pretty, encoding='UTF-8')

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
    buildXML(filename, roads, True)


if __name__ == "__main__":
    main()
