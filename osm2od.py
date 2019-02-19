import xml.etree.ElementTree
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
    e = xml.etree.ElementTree.parse('osmmaps/{}'.format(filename)).getroot()

    print("Reading file osmmaps/{}".format(filename))

    nodes = readNodes(e)
    roads = readRoads(e, nodes)

    print("Finished reading file, found {} nodes and {} roads.".format(len(nodes), len(roads)))
    

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('file', help="Input filename")

    args = parser.parse_args()

    if args.file:
        filename = args.file

    readOSM(filename)


if __name__ == "__main__":
    main()
