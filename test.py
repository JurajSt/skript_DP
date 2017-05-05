import ogr, os, math as m
import modul2
import definicie

lat = 49.03471077
lon = 20.32293131
elev = 746.018

lat_s = 49.03551077
lon_s = 20.32293131

hodnota = 4.44*m.pow(10,-4)
hodnota2 = m.pow(4.44*m.pow(10,-4),2)

xyz = modul2.fLatLonH_to_XYZ(lat, lon, elev)
x = xyz[0]
y = xyz[1]
z = xyz[2]
xyzs = modul2.fLatLonH_to_XYZ(lat_s, lon_s, elev)
xs = xyzs[0]
ys = xyzs[1]
zs = xyzs[2]
dx = x-xs
dy = y-ys

def fKruh(x, y, z):
    r = 65.0
    i = 0
    cesta_vystup_txt_kruh_xyz = "../data/vystup/kruznica_xyz.txt"
    txt_kruh_xyz = open(cesta_vystup_txt_kruh_xyz, "w")
    txt_kruh_xyz.write("ID X Y Z\n")
    cesta_vystup_txt_kruh_blh = "../data/vystup/kruznica_blh.txt"
    txt_kruh_blh = open(cesta_vystup_txt_kruh_blh, "w")
    txt_kruh_blh.write("ID L B H\n")
    idLH = 1
    idLD = 1
    idPH = 1
    idPD = 1
    blh = modul2.fXYZ_to_LatLonH(x,y,z)
    zapis_xyz = str(0) + " " + str(x).replace(".", ",") + " " + str(y).replace(".", ",") + " " + str(z).replace(".", ",") + "\n"
    zapis_blh = str(0) + " " + str(blh[1]).replace(".", ",") + " " + str(blh[0]).replace(".", ",") + " " + str(blh[2]).replace(".",",") + "\n"
    txt_kruh_xyz.write(zapis_xyz)
    txt_kruh_xyz.write(zapis_blh)
    LH = []
    LD = []
    PH = []
    PD = []
    while i<=r:
        yy = m.sqrt((r * r) - (i * i))
        ykLH = y + yy
        xkLH = x + i
        ykPH = y + yy
        xkPH = x - i
        ykPD = y - yy
        xkPD = x - i
        ykLD = y - yy
        xkLD = x + i
        blhLH = modul2.fXYZ_to_LatLonH(xkLH, ykLH, z)
        blhLD = modul2.fXYZ_to_LatLonH(xkLD, ykLD, z)
        blhPH = modul2.fXYZ_to_LatLonH(xkPH, ykPH, z)
        blhPD = modul2.fXYZ_to_LatLonH(xkPD, ykPD, z)
        zaznamLH = [idLH, xkLH, ykLH]
        zaznamLD = [idLD, xkLD, ykLD]
        zaznamPH = [idPH, xkPH, ykPH]
        zaznamPD = [idPD, xkPD, ykPD]
        LH.append(zaznamLH)
        LD.append(zaznamLD)
        PH.append(zaznamPH)
        PD.append(zaznamPD)
        zapis1 = str(idLH) + " " + str(xkLH).replace(".", ",") + " " + str(ykLH).replace(".", ",") \
                + " " + str(z).replace(".", ",") + "\n"
        zapis2 = str(idPH) + " " + str(xkPH).replace(".", ",") + " " + str(ykPH).replace(".", ",") \
                + " " + str(z).replace(".", ",") + "\n"
        zapis3 = str(idPD) + " " + str(xkPD).replace(".", ",") + " " + str(ykPD).replace(".", ",") \
                 + " " + str(z).replace(".", ",") + "\n"
        zapis4 = str(idLD) + " " + str(xkLD).replace(".", ",") + " " + str(ykLD).replace(".", ",") \
                 + " " + str(z).replace(".", ",") + "\n"
        zapis5 = str(idLH) + " " + str(blhLH[1]).replace(".", ",") + " " + str(blhLH[0]).replace(".", ",") \
                 + " " + str(blhLH[2]).replace(".", ",") + "\n"
        zapis6 = str(idLD) + " " + str(blhLD[1]).replace(".", ",") + " " + str(blhLD[0]).replace(".", ",") \
                 + " " + str(blhLD[2]).replace(".", ",") + "\n"
        zapis7 = str(idPH) + " " + str(blhPH[1]).replace(".", ",") + " " + str(blhPH[0]).replace(".", ",") \
                 + " " + str(blhPH[2]).replace(".", ",") + "\n"
        zapis8 = str(idPD) + " " + str(blhPD[1]).replace(".", ",") + " " + str(blhPD[0]).replace(".", ",") \
                 + " " + str(blhPD[2]).replace(".", ",") + "\n"
        txt_kruh_xyz.write(zapis1)
        txt_kruh_xyz.write(zapis2)
        txt_kruh_xyz.write(zapis3)
        txt_kruh_xyz.write(zapis4)
        txt_kruh_blh.write(zapis5)
        txt_kruh_blh.write(zapis6)
        txt_kruh_blh.write(zapis7)
        txt_kruh_blh.write(zapis8)
        idLH = idLH + 1
        idLD = idLD + 1
        idPH = idPH + 1
        idPD = idPD + 1
        i = i+0.01
    txt_kruh_xyz.close()
    txt_kruh_blh.close()
    return LH, LD, PH, PD

