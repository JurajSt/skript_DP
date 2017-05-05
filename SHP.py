import  ogr, os

def fKruhLinia(zoznam_suradnic, nazov, output):
    # Input data
    outSHPfn = output
    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(nazov, geom_type=ogr.wkbLineString)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    line = ogr.Geometry(ogr.wkbLineString)
    i = 1
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][1]
        Yd = zoznam_suradnic[i][2]
        # create point geometry
        line.AddPoint(Xd, Yd)
        i = i+1
    line.AddPoint(zoznam_suradnic[1][1], zoznam_suradnic[1][2]) # prvy bod uzaviera kruh
    outFeature.SetGeometry(line)
    outLayer.CreateFeature(outFeature)
    outFeature = None
    return 0

def fLineShp(zoznam_suradnic, nazov):
    # Input data
    fieldName = 'cas'
    fieldType = ogr.OFTInteger
    outSHPfn = "data/vystup/linia"
    # coordinate system
    #srs = osr.SpatialReference()   # hlada tabulku epsg kodov????
    #srs.ImportFromEPSG(4326)

    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(nazov, geom_type=ogr.wkbLineString)
    X = zoznam_suradnic[0][1]
    Y = zoznam_suradnic[0][2]
    i = 1
    # create a field
    idField = ogr.FieldDefn(fieldName, fieldType)
    outLayer.CreateField(idField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][1]
        Yd = zoznam_suradnic[i][2]
        fieldValue = zoznam_suradnic[i][0]
        # create point geometry
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(X, Y)
        line.AddPoint(Xd, Yd)
        i = i+1
        outFeature.SetGeometry(line)
        outFeature.SetField(fieldName, fieldValue)
        outLayer.CreateFeature(outFeature)
    outFeature = None
    return 0

def fKruhLinia(zoznam_suradnic, nazov, output):
    # Input data
    outSHPfn = output
    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(nazov, geom_type=ogr.wkbLineString)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    line = ogr.Geometry(ogr.wkbLineString)
    i = 1
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][1]
        Yd = zoznam_suradnic[i][2]
        # create point geometry
        line.AddPoint(Xd, Yd)
        i = i+1
    line.AddPoint(zoznam_suradnic[1][1], zoznam_suradnic[1][2]) # prvy bod uzaviera kruh
    outFeature.SetGeometry(line)
    outLayer.CreateFeature(outFeature)
    outFeature = None
    return 0

def fPointShp(zoznam_suradnic, nazov):
    # Input data
    fieldName = "azimut"
    fieldType = ogr.OFTInteger
    outSHPfn = "data/vystup/kruh/bod"
    n = "bod" + nazov
    # coordinate system
    #srs = osr.SpatialReference()   # hlada tabulku epsg kodov????
    #srs.ImportFromEPSG(4326)

    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(n, geom_type=ogr.wkbPoint)
    # create a field
    idField = ogr.FieldDefn(fieldName, fieldType)
    outLayer.CreateField(idField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    i = 0
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][0]
        Yd = zoznam_suradnic[i][1]
        fieldValue = zoznam_suradnic[i][2]
        # create point geometry
        point = ogr.Geometry(ogr.wkbPoint)

        point.AddPoint(Xd, Yd)
        i = i+1
        outFeature.SetGeometry(point)
        outFeature.SetField(fieldName, fieldValue)
        outLayer.CreateFeature(outFeature)
    outFeature = None
    return 0

def fLineClipShp(zoznam_suradnic, nazov):
    # Input data
    fieldName = 'cas'
    fieldName1 = 'druzica'
    fieldName2 = 'azimut'
    fieldType = ogr.OFTInteger
    fieldType1 = ogr.OFTString
    fieldType2 = ogr.OFTReal
    outSHPfn = "data/vystup/linia/orez"
    n = "orez_" + nazov
    # coordinate system
    #srs = osr.SpatialReference()   # hlada tabulku epsg kodov????
    #srs.ImportFromEPSG(4326)

    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(n, geom_type=ogr.wkbLineString)
    X = zoznam_suradnic[0][0]
    Y = zoznam_suradnic[0][1]
    i = 1
    # create a field
    idField = ogr.FieldDefn(fieldName, fieldType)
    idField1 = ogr.FieldDefn(fieldName1, fieldType1)
    idField2 = ogr.FieldDefn(fieldName2, fieldType2)
    outLayer.CreateField(idField)
    outLayer.CreateField(idField1)
    outLayer.CreateField(idField2)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][0]
        Yd = zoznam_suradnic[i][1]
        fieldValue = zoznam_suradnic[i][2]
        fieldValue1 = zoznam_suradnic[i][3]
        fieldValue2 = zoznam_suradnic[i][4]
        # create point geometry
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(X, Y)
        line.AddPoint(Xd, Yd)
        i = i+1
        outFeature.SetGeometry(line)
        outFeature.SetField(fieldName, fieldValue)
        outFeature.SetField(fieldName1, fieldValue1)
        outFeature.SetField(fieldName2, fieldValue2)
        outLayer.CreateFeature(outFeature)
    outFeature = None
    return 0