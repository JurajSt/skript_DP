import os
import sys
from simplekml import Kml
def fCreatePointKML(zoznam, n):
    kml = Kml()
    fol = kml.newfolder(name="kruh")
    At = zoznam[0][3]
    L = zoznam[0][0]
    B = zoznam[0][1]
    pnt = fol.newpoint(name="{0}".format(At), coords=[(L, B)])
    pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/paddle/ylw-blank-lv.png'
    pnt.style.iconstyle.scale = 0.5
    pnt.style.labelstyle.color = 'ffffffff'
    pnt.style.labelstyle.scale = 1
    for i in range(1,len(zoznam)):
        At = zoznam[i][3]
        if At == 0:
            continue
        L = zoznam[i][0]
        B = zoznam[i][1]
        pnt = fol.newpoint(name="{0}".format(At), coords=[(L, B)])
        pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/paddle/grn-blank-lv.png'
        pnt.style.iconstyle.scale = 0.4
        pnt.style.labelstyle.color = 'ffffffff'
        pnt.style.labelstyle.scale = 1
    kml.save("data/vystup/kruh/" + n + "kruh.kml")

def fCreateLineKML(zoznam, n):
    kml = Kml()
    fol = kml.newfolder(name="spojnica")
    At = zoznam[0][3]
    L = zoznam[0][0]
    B = zoznam[0][1]
    for i in range(1,len(zoznam)):
        At = zoznam[i][3]
        if At == 0:
            continue
        Ld = zoznam[i][0]
        Bd = zoznam[i][1]
        pnt = fol.newlinestring(name="{0}".format(At), coords=[(L, B), (Ld,Bd)])
    kml.save("data/vystup/kruh/" + n + "spojnica.kml")