# osm2opendrive
A tool for generating Apollo OpenDRIVE maps from OpenStreetMap data

## Apollo OpenDrive
Apollo is an architecture for autonomous driving developed by Baidu. They have made a modified version of the OpenDrive standard to suit their needs better. 

This project aims to create a tool for generating maps in the Apollo version of the OpenDrive specification using data from OpenStreetMap.

## Usage
```
usage: osm2od.py [-h] [-c CONFIG] [-z ZONE] [-p] file

positional arguments:
  file                  Input filename

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        Manually set lane numbers and widths based on road
                        names
  -z ZONE, --zone ZONE  UTM zone, example: -z 32V
  -p, --pretty          Prettify output
```

### Example config
Road name, Number of lanes, Lane width (meters)
```
Sem Sælands Vei, 1, 5
```
