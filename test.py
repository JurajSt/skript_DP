import ogr, os, math
import modul2
import definicie

lat = 49.03471077
lon = 20.32293131
elev = 746.018

lat_s = 49.03551077
lon_s = 20.32293131

hodnota = 4.44*math.pow(10,-4)
hodnota2 = math.pow(4.44*math.pow(10,-4),2)

la = lat_s - hodnota
lo = lon_s - hodnota2
print la, lo