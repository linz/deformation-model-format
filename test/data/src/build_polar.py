#!/usr/bin/python3
from LINZ.Geodetic.Ellipsoid import GRS80
import math
import sys
import csv

lons=[0.0,120.0,240.0,360.0]
lats=[-90.0,-89.5,-89.0]
en0=[1.3,2.5]
hfunc=lambda lon,lat: 3.0+(90+lat)*math.cos(math.radians(lon))*2.5

enu0=GRS80.enu_axes(0.0,-90.0)
evec=en0[0]*enu0[0]
nvec=en0[1]*enu0[1]
dvec=evec+nvec

def disp(lon,lat):
    enu=GRS80.enu_axes(lon,lat)
    de=enu[0].dot(dvec)
    dn=enu[1].dot(dvec)
    du=hfunc(lon,lat)
    return de,dn,du

points=[]
if len(sys.argv) > 1:
    csvr=csv.DictReader(open(sys.argv[1]))
    for r in csvr:
        points.append((float(r['lon']),float(r['lat'])))
else:
    for lat in lats: 
        for lon in lons:
            points.append((lon,lat))

print('lon,lat,de,dn,du')
for lon,lat in points:
    de,dn,du=disp(lon,lat)
    print('{0:.1f},{1:.1f},{2:.4f},{3:.4f},{4:.4f}'.format(lon,lat,de,dn,du))


