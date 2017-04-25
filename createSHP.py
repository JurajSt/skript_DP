import ogr, os
import modul2
lat = 49.03471077
lon = 20.32293131

# Input data
fieldName = 'line'
fieldType = ogr.OFTInteger
outSHPfn = 'test'
# Create the output shapefile
shpDriver = ogr.GetDriverByName("ESRI Shapefile")
if os.path.exists(outSHPfn):
    shpDriver.DeleteDataSource(outSHPfn)
outDataSource = shpDriver.CreateDataSource(outSHPfn)
outLayer = outDataSource.CreateLayer('clip', geom_type=ogr.wkbLineString)
# create point geometry
polygon = ogr.Geometry(type=ogr.wkbLineString)
x1 = modul2.fLatLonH_to_XYZ(49.0349077, 20.32333131, 746.018)
print x1
#polygon.AddPoint(20.32333131, 49.0349077)
polygon.AddPoint(x1[0], x1[1])
x1 = modul2.fLatLonH_to_XYZ(49.0349077, 20.32253131, 746.018)
print x1
#polygon.AddPoint(20.32253131, 49.0349077)
polygon.AddPoint(x1[0], x1[1])
x1 = modul2.fLatLonH_to_XYZ(49.03451384, 20.32253131, 746.018)
print x1
#polygon.AddPoint(20.32253131, 49.03451384)
polygon.AddPoint(x1[0], x1[1])
x1 = modul2.fLatLonH_to_XYZ(49.03451384, 20.32333131, 746.018)
print x1
#polygon.AddPoint(20.32333131, 49.03451384)
polygon.AddPoint(x1[0], x1[1])
x1 = modul2.fLatLonH_to_XYZ(49.0349077, 20.32333131, 746.018)
print x1
#polygon.AddPoint(20.32333131, 49.0349077)
polygon.AddPoint(x1[0], x1[1])
#create feature polygon:
feature = ogr.Feature(outLayer.GetLayerDefn() )
feature.SetGeometry(polygon)
outLayer.CreateFeature(feature)
feature = None