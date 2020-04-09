import gdal
import sys
import csv
import os.path
import argparse

parser=argparse.ArgumentParser(description='Dump deformation grid')
parser.add_argument('output_csv',help='Output CSV file')
parser.add_argument('deformation_geotiff',nargs='+',help='Deformation GeoTIFF file')
args=parser.parse_args()

if not args.output_csv.endswith('.csv'):
    raise RuntimeError('Output must be a .csv file')

with open(args.output_csv,'w') as csvh:
    csvw=csv.writer(csvh)
    csvw.writerow('filename subgrid name parent nchild dtype utype dx dy wkt'.split())
    for gtiff in args.deformation_geotiff:
        tfile=os.path.basename(gtiff)
        ds=gdal.Open(gtiff)
        subds=ds.GetSubDatasets()
        if len(subds) == 0:
            subds=[(gtiff,None)]
        for ds,ignore in subds:
            # print("Loading {0}".format(ds))
            grid=gdal.Open(ds)
            ds=ds.replace(gtiff,'')
            md=grid.GetMetadata()
            name=md.get('grid_name','')
            parent=md.get('parent_grid_name','')
            nchild=md.get('number_of_nested_grids','')
            dtype=md.get('DISPLACEMENT_TYPE','').lower()
            utype=md.get('UNCERTAINTY_TYPE','').lower()
            gt=grid.GetGeoTransform()
            x0=round(gt[0],10)
            y0=round(gt[3],10)
            dx=round(gt[1],10)
            dy=round(gt[5],10)
            nx=grid.RasterXSize
            ny=grid.RasterYSize
            x0 += dx/2
            y0 += dy/2
            x1 = x0+dx*nx
            y1 = y0+dy*ny
            shape='POLYGON(({0} {1},{0} {3},{2} {3},{2} {1},{0} {1}))'.format(
                x0,y0,x1,y1)
            csvw.writerow((tfile,ds,name,parent,nchild,dtype,utype,str(dx),str(dy),shape))






