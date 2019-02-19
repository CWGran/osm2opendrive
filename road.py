class Road:
    
    id = ""
    nodes = []

    def __init__(self, id):
        self.id = id

class Node:
    
    id = ""
    lat = ""
    lng = ""

    def __init__(self, id, lat, lng):
        self.id = id
        self.lat = lat
        self.lng = lng
