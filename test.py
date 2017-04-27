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
    cesta_vystup_txt_kruh = "../data/vystup/kruznica.txt"
    txt_poloha_d = open(cesta_vystup_txt_kruh, "w")
    txt_poloha_d.write("ID B L H\n")
    idLH = 1
    idLD = 1
    idPH = 1
    idPD = 1
    zapis = str(0) + " " + str(x).replace(".", ",") + " " + str(y).replace(".", ",") + " " + str(z).replace(".", ",") + "\n"
    txt_poloha_d.write(zapis)
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
        #zapis = str(id) + " " + str(blhk[0]).replace(".", ",") + " " + str(blhk[1]).replace(".", ",") + " " + str(blhk[2]).replace(".", ",") + "\n"
        txt_poloha_d.write(zapis1)
        txt_poloha_d.write(zapis2)
        txt_poloha_d.write(zapis3)
        txt_poloha_d.write(zapis4)
        idLH = idLH + 1
        idLD = idLD + 1
        idPH = idPH + 1
        idPD = idPD + 1
        i = i+0.01
    txt_poloha_d.close()
    return LH, LD, PH, PD

