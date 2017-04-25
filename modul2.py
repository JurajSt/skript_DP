import math
import os
import sys
import ogr
import osgeo.osr as osr
#-------------------------------------------------------------------------------
#definuje nazov podla vstupnych suborov
def fnazov(cesta):
    pocetZnakov = 0
    for znak in cesta:
        if "/" in znak:
            pocetZnakov = pocetZnakov+ 1
    nazov = cesta.split("/")[pocetZnakov].split(".")[0]
    return nazov
#-------------------------------------------------------------------------------
#prevod XYZ do geografickych souradnic, pokud se vezmou souradnice z TRO, presnost
#zpetne trasnformace do XYZ je dle testovani do 1 mm ve vsech slozkach souradnic
#postup prevzat z Hoffman-Wellenhof, provadi se dve iterace
def fXYZ_to_LatLonH(X,Y,Z):
    X=float(X)
    Y=float(Y)
    Z=float(Z)

    a=6378137.0
    b=6356752.3142
    e_2=((a*a)-(b*b))/(a*a)
    e=math.sqrt(e_2)

    lon = math.atan(Y/X)
    lon_degrees = math.degrees(lon)

    p = math.sqrt(X*X+Y*Y)

    lat0=math.atan((Z/p)*(1/(1-e_2)))
    cos2lat0 = (1+math.cos(2*lat0))/2
    sin2lat0 = (1-math.cos(2*lat0))/2
    N0 = (a*a)/math.sqrt(a*a*cos2lat0 + b*b*sin2lat0)
    h0=(p/math.cos(lat0))-N0
    lat1=math.atan((Z/p)*(1/(1-(e_2)*(N0/(N0+h0)))))

    cos2lat1 = (1+math.cos(2*lat1))/2
    sin2lat1 = (1-math.cos(2*lat1))/2
    N1 = (a*a)/math.sqrt(a*a*cos2lat1 + b*b*sin2lat1)
    h1=(p/math.cos(lat1))-N1
    lat2=math.atan((Z/p)*(1/(1-(e_2)*(N1/(N1+h1)))))
    lat2_degrees = math.degrees(lat2)

    cos2lat2 = (1+math.cos(2*lat2))/2
    sin2lat2 = (1-math.cos(2*lat2))/2
    N2 = (a*a)/math.sqrt(a*a*cos2lat2 + b*b*sin2lat2)
    h2=(p/math.cos(lat2))-N2
    lat3=math.atan((Z/p)*(1/(1-(e_2)*(N2/(N2+h2)))))
    lat3_degrees = math.degrees(lat3)

    if X < 0:
        lon_degrees = lon_degrees + 180

    return lat3_degrees, lon_degrees, h2

def fLatLonH_to_XYZ(B,L,H):

    a=6378137.0
    b=6356752.3142
    e_2=((a*a)-(b*b))/(a*a)
    B = math.radians(float(B))
    L = math.radians(float(L))

    N = a/math.sqrt(1-e_2*math.pow(math.sin(B),2))

    X = (N+H)*math.cos(B)*math.cos(L)
    Y = (N+H)*math.cos(B)*math.sin(L)
    Z = (N*(1-e_2)+H)*math.sin(B)

    return X, Y, Z

# zistuje uhol medzi priamkou ktora je tvorena suradnicami stanice a
# suradnicami stanice prepocitane z BL (bez elipsoidickej vysky),
# druha priamka je tvorena suradnicami stanice a druzice
def fUhol_priamok(Bod_X,Bod_Y,Bod_Z, X, Y, Z, ):

    dx = X-Bod_X
    dy = Y-Bod_Y
    dz = Z-Bod_Z
    dlzka = math.sqrt(math.pow(dx,2) + math.pow(dy,2) + math.pow(dz,2))

    bodBLH = fXYZ_to_LatLonH(Bod_X, Bod_Y, Bod_Z)
    Bod_B = bodBLH[0]
    Bod_L = bodBLH[1]
    bod_elipsoid = fLatLonH_to_XYZ(Bod_B, Bod_L, 0)

    dx_bod = bod_elipsoid[0] - Bod_X
    dy_bod = bod_elipsoid[1] - Bod_Y
    dz_bod = bod_elipsoid[2] - Bod_Z
    dlzka_bod = math.sqrt(math.pow(dx_bod,2) + math.pow(dy_bod,2) + math.pow(dz_bod,2))

    uhol_rad = -1*((dx*dx_bod + dy*dy_bod + dz*dz_bod)/(dlzka*dlzka_bod))
    uhol_stupne = math.degrees(uhol_rad)

    return uhol_stupne

def fvypocet_poloha(g, Dt):
    # konstanty
    sprava = []
##    print g
    omega = 0.000072921151467           # uhlova rychlost rotacie zeme
    mi = 398600500000000.0              # geocentricka gravitacna konstanta
    pi = math.pi

    for riadok in g:
        riadok0 = riadok.replace("D-","De")             # replace - nahradza v rinexe, prikald: .515379774666D+04 => 0.515379774666e+04
        ##riadok1 = riadok0.replace("D","e")
        riadok2 = riadok0.replace("-"," -")
        riadok3 = riadok2.replace("De","e-")
        riadok4 = riadok3.replace("D","e")
        riadok5 = riadok4.replace("Em","E-")
        sprava.append(riadok5)
        #print sprava
##    Toe = float(sprava[3].split()[0])
##    print Toe
##    dodatok_cas_nav = (Toe%60)%60
##    dodatok_cas_nav1 = 60 - dodatok_cas_nav
##    if dodatok_cas_nav1 < 60.0:
##        Dt = Dt + dodatok_cas_nav1
##        print Dt, dodatok_cas_nav
    n0 = math.sqrt((mi / math.pow(math.pow(float(sprava[2].split()[3]),2),3)))
    n = n0 + float(sprava[1].split()[2])
    M = float(sprava[1].split()[3]) + n * Dt
##    print float(sprava[1].split()[3])
##    E = M
    a = M
    b = 1

    while round(a,5) != round(b,5):
        E0= a - ((a-float(sprava[2].split()[1]) * math.sin(a) - M)/(1 - float(sprava[2].split()[1]) * math.cos(a)))
        E1 = E0 - ((E0-float(sprava[2].split()[1]) * math.sin(E0) - M)/(1 - float(sprava[2].split()[1]) * math.cos(E0)))
        del(a,b)                # mazem z dovodu prepisovania premennej
        a = E1
        b = E0
##        print a, b

    cosv = (math.cos(E1)-float(sprava[2].split()[1])) / (1-float(sprava[2].split()[1]) * math.cos(E1))
    if E1<0:
        E1 = (2*pi)+E1
    if E1>pi:
        v =(2*pi)-math.acos(cosv)
    else:
        v = math.acos(cosv)
    w = float(sprava[4].split()[2])

    cos2x = math.pow(math.cos(v+w),2) - math.pow(math.sin(v+w),2)
    sin2x = 2 * math.sin(v+w) * math.cos(v+w)

    du = float(sprava[2].split()[0]) * cos2x + float(sprava[2].split()[2]) * sin2x
    u = v + w + du

    dr = float(sprava[4].split()[1]) * cos2x + float(sprava[1].split()[1]) * sin2x
    r = math.pow(float(sprava[2].split()[3]),2) * (1 - float(sprava[2].split()[1]) * math.cos(E1)) + dr

    x_p = r * math.cos(u)
    y_p = r * math.sin(u)
    r_kontrola = math.sqrt(math.pow(x_p,2) + math.pow(y_p,2))

##    print "r", r,"r_konnt", r_kontrola
##    if str(r) == str(r_kontrola):                   # taktiez musim dat str() inak padmienka nefunguje
##            print "vypocet je zatial v poriadku"    # podminka kontroluje spravnost vypoctu
##            print r, r_kontrola
##    elif str(r) != str(r_kontrola):                 # taktiez musim dat str inak padmienka nefunguje
##        print "vo vypocte je chyba"                 # bez str() urci dne rovnake hodnoty ako nerovne
##        print r, r_kontrola                         # a skript skonci
##        os.system("pause")
##        sys.exit()

    di = float(sprava[3].split()[1]) * cos2x + float(sprava[3].split()[3]) * sin2x
    i = float(sprava[4].split()[0]) + float(sprava[5].split()[0]) * Dt + di
    O = float(sprava[3].split()[2]) + (float(sprava[4].split()[3]) - omega) * Dt - omega * float(sprava[3].split()[0])

    X = x_p * math.cos(O) - y_p * math.sin(O) * math.cos(i)
    Y = x_p * math.sin(O) + y_p * math.cos(O) * math.cos(i)
    Z = y_p * math.sin(i)

##    pole_XYZ.append(sprava[0][0:2].replace(" ",""))
##    pole_XYZ.append(cas_spravy)
##    pole_XYZ.append(X)
##    pole_XYZ.append(Y)
##    pole_XYZ.append(Z)

    return X, Y, Z

def fData_EPH(telo_eph, cd_nav, cas_sekundy_nav):
    for data_eph in telo_eph:
        if "EOF" in data_eph:
            konec = "konec zaznamu"
            return "konec zaznamu"
        if "*" in data_eph:
            data_eph1 = data_eph.split()
##            print data_eph1
            rok_eph = int(data_eph1[1])
##            print data_eph1[1]
            mes_eph = int(data_eph1[2])
            den_eph = int(data_eph1[3])
            hod_eph = int(data_eph1[4])
            min_eph = int(data_eph1[5])
            sek_eph = float(data_eph1[6])
            cas_sekundy_eph = (hod_eph*60*60) + (min_eph*60) + sek_eph

            continue
        data_eph1 = data_eph.split()
        cd_eph = data_eph1[0]
##        print cd_eph
        if "PR" in cd_eph:
            continue
        if cas_sekundy_eph == cas_sekundy_nav:
            if cd_eph == cd_nav:
##                print cas_sekundy_eph
                return data_eph1, cas_sekundy_eph, rok_eph, mes_eph, den_eph

def fEph_zapis(telo_eph, nazov, nazov_eph, koncovka_txt):
    # txt - vystup pre efemeridy
    cesta_vystup_eph_txt = "../data/vystup/" + nazov + nazov_eph + koncovka_txt
    txt_eph = open(cesta_vystup_eph_txt, "w")
    txt_eph.write("Cislo_druzice cas_eph X_ Y_ Z_\n")
    for data_eph in telo_eph:
        if "EOF" in data_eph:
            return "konec zaznamu"
        if "*" in data_eph:
            data_eph1 = data_eph.split()
            hod_eph = int(data_eph1[4])
            min_eph = int(data_eph1[5])
            sek_eph = float(data_eph1[6])
            cas_eph = (hod_eph * 60 * 60) + (min_eph * 60) + sek_eph
            continue
        data_eph1 = data_eph.split()
        cd_eph = data_eph1[0]
        X_eph = float(data_eph1[1]) * 1000
        Y_eph = float(data_eph1[2]) * 1000
        Z_eph = float(data_eph1[3]) * 1000
        if "PR" in cd_eph:
            continue
        zapis_eph_txt = cd_eph + " " + str(cas_eph) + " " + str(X_eph) + " " + str(Y_eph) + " " + str(
            Z_eph) + "\n"  # zapis vyslednych hodnot txt
        txt_eph.write(zapis_eph_txt)
    txt_eph.close()

def fObservacie(telo, pocet_kan):
    ttt = telo[1].strip("\n")
    if "COMMENT" in ttt:     # v kazdej hodine su vlozene komenty krore umazavam
        del(telo[:1])
        while "COMMENT" in telo[0]:
            del (telo[0])
    riadok1 = telo[0].split()
    rok_obs = int(riadok1[0])
    mes_obs = int(riadok1[1])
    den_obs = int(riadok1[2])
    hod_obs = int(riadok1[3])
    min_obs = int(riadok1[4])
    sek_obs = float(riadok1[5])
    cas_obs = ((hod_obs*60)*60) + (min_obs*60) + sek_obs
    zaznam_cas_obs = []
    zaznam_cas_obs.append(rok_obs)
    zaznam_cas_obs.append(mes_obs)
    zaznam_cas_obs.append(den_obs)
    zaznam_cas_obs.append(cas_obs)
    pocet_druzic = int(riadok1[len(riadok1) - 1][:2])  # pocet druzic
    pocet_cisel_druzic = (len(riadok1[len(riadok1) - 1][2:]) / 3)  # pocet cisiel druzic v prvom riadku
    podiel1 = pocet_druzic % pocet_cisel_druzic # delim int - cislo je odstrihnute od ciarky, +1 pre indexaciu v cykle
    podiel2 = pocet_druzic / pocet_cisel_druzic  # delim int - cislo je odstrihnute od ciarky, +1 pre indexaciu v cykle
    if podiel1>0:
        podiel = podiel2+1
    else:
        podiel = podiel2
    cisla_druzic = riadok1[len(riadok1) - 1][2:]
    for data in telo[1:podiel]:
        riadok2 = data.split()
        cisla_druzic = cisla_druzic + riadok2[0]  # mam vsetky cisla druzic v jedonom riadku
    poc_riadkov_obs = (pocet_kan / 5) + 1
    zaznam = ""
    zaznam_cas_obs.append(cisla_druzic)
    observacie = []
    i = 0
    for data2 in telo[podiel:podiel + poc_riadkov_obs * pocet_druzic]:
        riadok3 = data2.strip("\n")
        while len(riadok3)<80:
            riadok3 = riadok3 + " "
        zaznam = zaznam + riadok3
        i = i + 1
        if i == poc_riadkov_obs:
            observacie.append(zaznam)
            zaznam = ""
            i = 0
    del telo[:podiel + poc_riadkov_obs * pocet_druzic]  # umazem spracovanu cas dat
    return zaznam_cas_obs, observacie

# vyhldanie v ascii suboroch vysku elipsoidu nad geoidom
def fhel_to_geoid(B, L):

    B = float(B)
    L = float(L)

    n = -90
    s = -45
    for i in range(1,5):
        if int(B) in range(n,s):
##           print n
           break
        n = n+45
        s = s+45

    w = -180
    e = -135
    for i in range(1,9):
        if int(L) in range(w,e):
##              print w
              break
        w = w+45
        e = e+45

    if n >= 0 and w >= 0:
        nazov ="n"+str(abs(n))+"e"+str(abs(w))
        nazov = nazov.replace("n0","n00")
        nazov = nazov.replace("e0","e00")
    elif n >= 0 and w < 0:
        nazov ="n"+str(abs(n))+"w"+str(abs(w))
        nazov = nazov.replace("n0","n00")
    elif n < 0 and w < 0:
        nazov ="s"+str(abs(n))+"w"+str(abs(w))
    elif n < 0 and w >= 0:
        nazov ="s"+str(abs(n))+"e"+str(abs(w))
        nazov = nazov.replace("e0","e00")
##    print nazov

    cesta_raster = os.path.join("../data/vstup/geoid_2008/ascii/"+nazov+".asc")
    if not os.path.exists(cesta_raster):
        print " ascii subor pre geoid neexistuje, navratova hodnota = 0"
        vyska = 0
        return vyska
##        sys.exit()

    vstup_raster = open(cesta_raster, "r")
    subor_raster = vstup_raster.readlines()
    vstup_raster.close()

    hlavicka = subor_raster[0:6]
    raster = subor_raster[6:]
    ##print hlavicka
    #print len(raster)

    stlpec = int(hlavicka[0].split()[1])
    riadok = int(hlavicka[1].split()[1])
    x_zaciatok = int(hlavicka[2].split()[1])
    y_zaciatok = int(hlavicka[3].split()[1])
    velkost_b = float(hlavicka[4].split()[1].replace(",","."))
    ##print stlpec, riadok, x_zaciatok, y_zaciatok, velkost_b

    vypocet_stlpca = (int(round((L-x_zaciatok)/velkost_b,0)))-1
    vypocet_riadka = (riadok-abs(int(round((B-y_zaciatok)/velkost_b,0))))-1
#    print vypocet_riadka, vypocet_stlpca
    for riadok_r in range(len(raster)):
        if riadok_r == vypocet_riadka:
            hodnota = raster[vypocet_riadka].split()[vypocet_stlpca]
            #print hodnota
            break
        #if riadok == vypocet_riadka:
            #hodnota = raster[vypocet_riadka-1].split()[vypocet_stlpca]
            #break
    vyska = float(hodnota.replace(",","."))
    del(subor_raster, hlavicka, raster)
    return vyska

def fCalculateAzimuth(xf,yf, xl, yl):
    dX = xl - xf
    dY = yl - yf
    PI = math.pi
    Azimuth = 0  # Default case, dX = 0 and dY >= 0
    if dX > 0:
        Azimuth = 90 - math.atan(dY / dX) * 180 / PI
    elif dX < 0:
        Azimuth = 270 - math.atan(dY / dX) * 180 / PI
    elif dY < 0:
        Azimuth = 180

    return Azimuth


def fLineShp(zoznam_suradnic, nazov):
    # Input data
    fieldName = 'cas'
    fieldType = ogr.OFTInteger
    outSHPfn = 'shp'
    # coordinate system
    #srs = osr.SpatialReference()   # hlada tabulku epsg kodov????
    #srs.ImportFromEPSG(4326)

    # Create the output shapefile
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outSHPfn):
        shpDriver.DeleteDataSource(outSHPfn)
    outDataSource = shpDriver.CreateDataSource(outSHPfn)
    outLayer = outDataSource.CreateLayer(nazov, geom_type=ogr.wkbLineString)
    X = zoznam_suradnic[0][2]
    Y = zoznam_suradnic[0][1]
    i = 1
    # create a field
    idField = ogr.FieldDefn(fieldName, fieldType)
    outLayer.CreateField(idField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    while i < len(zoznam_suradnic):
        Xd = zoznam_suradnic[i][2]
        Yd = zoznam_suradnic[i][1]
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

def fSkratenie(x, y, xd, yd):
    dd = 20
    dx = x - xd
    dy = y - yd
    d = math.sqrt((dx * dx) + (dy * dy))
    a = dy / d
    b = dx / d
    dyy = a * dd
    dxx = b * dd
    xxd = x + dxx
    yyd = y + dyy
    return xxd, yyd