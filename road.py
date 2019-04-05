class Road:
    
    def __init__(self, id):
        self.id = id
        self.nodes = []

class Node:
    
    def __init__(self, id, lat, lng):
        self.id = id
        self.lat = float(lat)
        self.lng = float(lng)
