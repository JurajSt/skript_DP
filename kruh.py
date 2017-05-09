import math as m
import SHP, modul2

def fkruznica(lat, lon, n):         # n je nazov stanice
    nazov = "kruznica_b_" + n
    cesta_vystup_txt_kruh_blh = "data/vystup/kruh/" + nazov + ".txt"
    txt_kruh_blh = open(cesta_vystup_txt_kruh_blh, "w")
    txt_kruh_blh.write("Id L B A At\n")
    r = 65                  # polomer kruhu x y
    rozdiel = 0.0008        # polomer kruhu lat lon
    i=0
    id = 0
    zoznam = []
    vertex = []
    while i <= 90:
        alfa = i            # alfa v stupnoch
        alfa_rad = m.radians(alfa)
        a = m.sin(alfa_rad)*r
        b = m.cos(alfa_rad)*r
        xx = (b*rozdiel)/r
        yy = (a*rozdiel)/r
        xy = [xx,yy]
        zoznam.append(xy)
        i = i+1
    zoznambl = []
    for j in range(len(zoznam)):
        x = lon + zoznam[j][0]
        y = lat + zoznam[j][1]
        bl = [x,y]
        v = [id, x, y]
        zoznambl.append(bl)
        vertex.append(v)
        id = id + 1
        #zapis1 = str(id) + " " + str(x).replace(".", ",") + " " + str(y).replace(".", ",") + "\n"
        #txt_kruh_blh.write(zapis1)
    zoznam.reverse()
    for j in range(len(zoznam)):
        x = lon - zoznam[j][0]
        y = lat + zoznam[j][1]
        bl = [x,y]
        v = [id, x, y]
        zoznambl.append(bl)
        vertex.append(v)
        id = id + 1
        #zapis1 = str(id) + " " + str(x).replace(".", ",") + " " + str(y).replace(".", ",") + "\n"
        #txt_kruh_blh.write(zapis1)
    zoznam.reverse()
    for j in range(len(zoznam)):
        x = lon - zoznam[j][0]
        y = lat - zoznam[j][1]
        bl = [x,y]
        v = [id, x, y]
        zoznambl.append(bl)
        vertex.append(v)
        id = id + 1
        #zapis1 = str(id) + " " + str(x).replace(".", ",") + " " + str(y).replace(".", ",") + "\n"
        #txt_kruh_blh.write(zapis1)
    zoznam.reverse()
    for j in range(len(zoznam)):
        x = lon + zoznam[j][0]
        y = lat - zoznam[j][1]
        bl = [x,y]
        v = [id, x, y]
        zoznambl.append(bl)
        vertex.append(v)
        id = id + 1
        #txt_kruh_blh.write(zapis1)

    j=0
    while j < 360:
        pocet = zoznambl.count(zoznambl[j])
        if pocet > 1:
            zoznambl.remove(zoznambl[j])
        j = j + 1
    for i in range(len(zoznambl)):
        azimut = modul2.fCalculateAzimuth(lon, lat, zoznambl[i][0], zoznambl[i][1])
        zoznambl[i].append(azimut)
        if azimut == 0:
            At = "S"
            zoznambl[i].append(At)
        elif azimut == 90:
            At = "V"
            zoznambl[i].append(At)
        elif azimut == 180:
            At = "J"
            zoznambl[i].append(At)
        elif azimut == 270:
            At = "Z"
            zoznambl[i].append(At)
        elif round(azimut,0)%10 == 0:
            At = int(round(azimut,0))
            zoznambl[i].append(At)
        else:
            At = 0
            zoznambl[i].append(At)
        zapis1 = str(id) + " " + str(zoznambl[i][0]).replace(".", ",") + " " + str(zoznambl[i][1]).replace(".", ",") + " " +\
                 str(azimut).replace(".", ",") + " " + str(At) + "\n"
        txt_kruh_blh.write(zapis1)
    zoznambl.insert(0, [lon, lat, 0, n])
    txt_kruh_blh.write("0" + " " + str(lon).replace(".", ",") + " " + str(lat).replace(".", ",") + " 0 " + n +"\n" )
    cestaSHP = "data/vystup/kruh"
    linia = SHP.fKruhLinia(vertex, nazov, cestaSHP)
    #bodKruh = SHP.fPointShp(zoznambl, nazov)
    txt_kruh_blh.close()
    return zoznambl