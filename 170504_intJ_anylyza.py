import sys
import os
import modul2, SHP, kruh, priesecnik
import math
import numpy as np
pi = math.pi

# vstupy
# presne efemeridy
cesta_eph = os.path.join("../data/vstup/cod19212.eph")  # cof19196.eph
if not os.path.exists(cesta_eph):
    print "subor eph neexistuje"
    sys.exit()
## navigacna sprava
cesta_nav = os.path.join("../data/vstup/ganp3060.16n")  # gan60.16n gope0480.14n
if not os.path.exists(cesta_nav):
    print "subor nav neexistuje"
    sys.exit()
# observacne data
cesta_obs = os.path.join("../data/vstup/ganp3060.16o")
if not os.path.exists(cesta_obs):
    print "subor neexistuje"
    sys.exit()
vstup_eph = open(cesta_eph, "r")
vstup_nav = open(cesta_nav, "r")
vstup_obs = open(cesta_obs, "r")
subor_eph = vstup_eph.readlines()
subor_nav = vstup_nav.readlines()
subor_obs = vstup_obs.readlines()
vstup_eph.close()
vstup_nav.close()
vstup_obs.close()
# ###########################################################
#  vystupy
zoznam_farieb = 'ff4672eb', 'ffda6e6f', 'ff9b26d7', 'ff941645', 'ffff220f', 'ffe09ea7', 'ff0a0a11', 'ff707236', 'ffc6cc23', 'ff2c8424', 'fff10b5a', 'ffcee25d', 'ff5538f2', 'fff708ae', 'ff5d18d2', 'fff50f2e', 'ffef66c1', 'ff464eec', 'ff354725', 'ffb1b329', 'ff643f30', 'ffe5b0fd', 'ffe512a6', 'ffaa0d2c', 'ffdd618a', 'ff57aa6a', 'ffa866da', 'ff430189', 'ff50757b', 'ff93568b', 'ff504621', 'ff6374b9', 'ff9e1e08', 'ffc52878', 'ff2a1c75', 'ff603267', 'ffa1ba0f', 'ff24a888', 'ffa0c757', 'ff413393', 'ffcadde6', 'ffb1bc4c', 'ffc5962e', 'ff9e9e60', 'ff3fcdc7', 'ff95bb90', 'ffbe26b5', 'ff81d106', 'ff1b8d9b', 'fffcde55', 'ffebbf0f', 'ff5caaa7', 'ff07bdd1', 'ffa36485', 'ff2ffc4d', 'ffec06e7', 'ff8eec68', 'ff46a9de', 'ff844149', 'fff323e2', 'ff433804', 'ff6caaec', 'ffc9392a', 'fff3543e', 'ff2d2755', 'ff42b097', 'ffe89bef', 'ff550801', 'ffddefdb', 'ff86d370', 'ff8170a3', 'ff2f84ea', 'ff04d250', 'ff433845', 'ffb22ca7', 'ffea79b4', 'ff865deb', 'ff8fba95', 'ff83ecdc', 'ff84726a', 'ff635c0b', 'ffb6cb0f', 'ff782274', 'ffe984b1', 'ff90a410', 'ffd72f38', 'ff17607b', 'ffebd493', 'ff028e9b', 'ff048cdd', 'ff1b5be7', 'ff144a1d', 'ff6f60d4', 'ff5079df', 'ff7987f7', 'ff31076d', 'ff7a9827', 'ffa83f21', 'ff3c63b4', 'ff9aca49', 'ff101008'
koncovka_xls = ".xls"
koncovka_txt = ".txt"
nazov = modul2.fnazov(cesta_nav)  # funkcia vracia nazov vstupneho suboru
nazov_eph = modul2.fnazov(cesta_eph)
# xls - porovnanie eph a vypocitane suradnice druzice
cesta_vystup_xls_porovnanie = "../data/vystup/" + nazov + nazov_eph + "_porovnanie" + koncovka_xls
xls_porovnanie = open(cesta_vystup_xls_porovnanie, "w")
# txt - poloha druzic
cesta_vystup_txt_poloha_d = "../data/vystup/" + nazov + koncovka_txt
txt_poloha_d = open(cesta_vystup_txt_poloha_d, "w")
# xls - vytup s observacnymi datami priradenych druzici + suradnice XYZ, BLH
cesta_vystup_xls_obs = "../data/vystup/" + nazov + koncovka_xls
xls_obs = open(cesta_vystup_xls_obs, "w")

# zapis haviciek
cesta_vystup_txt_orez = "../data/vystup/orez_" + nazov + koncovka_txt
txt_orez = open(cesta_vystup_txt_orez, "w")
xls_obs.write("CD\tcas\tsin_elevUhol\tstupne\tSRN1\tSRN_lin1\tSRN2\tSRN_lin2\tXd\tYd\tZd\tBd\tLd\tHd\tazimut\td\n")
xls_porovnanie.write(
    "CD_eph\tCD_nav\tcas_eph\tcas_nav\tDt\tXeph\tYeph\tZeph\tXvyp\tYvyp\tZvyp\tXroz\tYroz\tZroz\tVektor\n")
txt_poloha_d.write("Cislo_druzice cas_eph cas_nav Dt X_vyp Y_vyp Z_vyp\n")
txt_orez.write("ID X Y\n")
# ###########################################################
# nacitanie udajov z eph
i = 0
while "*" is not subor_eph[i][0]:
    i = i + 1
telo_eph = subor_eph[i:]
# ##########################################################
# nacitanie udajov z hlavicky observacii
hlavicka_obs = []
for data in subor_obs:
    hlavicka_obs.append(data)
    if "END" in data:
        ##        print data
        break

    ts1 = ""
    TappPosXYZ = "APPROX POSITION XYZ"
    Tant = "ANTENNA: DELTA H/E/N"
    Tint = "INTERVAL"
    TObs = "# / TYPES OF OBSERV"
    TfObs = "TIME OF FIRST OBS"


    ##print hlavicka
    for data in hlavicka_obs:
        if "APPROX POSITION XYZ" in data:
            appPosXYZ = data[0:(len(data) - (len(TappPosXYZ) + 1))]
        if "ANTENNA: DELTA H/E/N" in data:
            antHEN = data[0:(len(data) - (len(Tant) + 1))]
        if "INTERVAL" in data:
            interval = float(data[0:(len(data) - (len(Tint) + 1))])
        if "# / TYPES OF OBSERV" in data:
            typObs = data[0:(len(data) - (len(TObs) + 1))]
            ts = typObs + " "
            ts1 = ts1 + ts
        if "TIME" in data:
            firstObs = data[0:(len(data) - (len(TfObs) + 1))]
        if "latitude" in data:
            Latitude = data.split()
            Lat = Latitude[0]
        if "longitude" in data:
            Longitude = data.split()
            Lon = Longitude[0]
        if "elevation" in data:
            Elevation = data.split()
            Elev = Elevation[0]

first_observ = firstObs.split()
appPos_XYZ = appPosXYZ.split()
ant_HEN = antHEN.split()
typ_observ = ts1.split()
pocet_kan = int(typ_observ[0])
typ_observ.remove(typ_observ[0])
# print "pocet kanalov:", pocet_kan
# print "suradnice XYZ:", appPos_XYZ
# print "pocet a signaly:", typ_observ
# print "prva observacia:", first_observ
# print "interval observacii", interval
# print "vysky anteny", ant_HEN
pilier = 2.3               # vyska piliera
B = float(Lat)
L = float(Lon)
H = float(Elev)
XYZ = modul2.fLatLonH_to_XYZ(0, L, H)
X = XYZ[0]
Y = XYZ[1]
Z = XYZ[2]
#H_geoid = H - modul2.fhel_to_geoid(B, L)
dr = "PG23"
## ##########################################################
## nacitanie udajov z navigacnej spravy
hlavicka_nav = []
for data in subor_nav:
    hlavicka_nav.append(data)
    if "END" in data:
        break
telo_nav = subor_nav[len(hlavicka_nav):]
# print telo_nav
hodnoty = []
q = 0
while q < len(telo_nav):  # rozdelenie tela na jednotlive spravy
    hodnoty.append(telo_nav[q:q + 8])
    q = q + 8
# print hodnoty
platnost_spravy = 7200  # 7200 sekund = 2 hodiny
k = 0

vektor_max = 0
vektor_min = 9999999
kontrola = modul2.fEph_zapis(telo_eph, nazov, nazov_eph, koncovka_txt)
data_all = []
while k < len(hodnoty):
    sprava = hodnoty[k]  # nacitanie jednej spravy
    #     print sprava
    if len(sprava) < 8:
        k = k + 1
        continue

    cd_nav = sprava[0].split()[0]  # cislo druzice
    if len(cd_nav) == 1:
        cd_nav = "PG0" + cd_nav
    else:
        cd_nav = "PG" + cd_nav
        # if cd_nav != dr:
        # k = k + 1
        # continue

    rok_nav = int(sprava[0].split()[1])  # urcenie datumu spravy
    mes_nav = int(sprava[0].split()[2])
    den_nav = int(sprava[0].split()[3])  # prepocet casu z navig. spravy na sekundy
    hod_nav = int(sprava[0].split()[4])
    min_nav = int(sprava[0].split()[5])
    sek_nav = float(sprava[0].replace("-", " -").split()[6])
    cas_nav = (hod_nav * 60 * 60) + (min_nav * 60) + sek_nav

    dodatok = (60 - sek_nav)
    Dt = 0
    if dodatok < 60:  # niektore spravy su nemaju celi cas napr: 11:59:44... do vypoctu treba zahrnut
        Dt = Dt + dodatok
    else:
        dodatok = 0
    while Dt < (platnost_spravy + dodatok):
        data = []
        poloha = modul2.fvypocet_poloha(sprava, Dt)  # vypocet polohy druzice v zadanom case Dt
        # EPH = modul2.fData_EPH(telo_eph, cd_nav, (cas_nav + Dt))
        # if EPH == "konec zaznamu":
        #    break
        # rok_eph = int(str(EPH[2])[2:])
        # mes_eph = EPH[3]
        # den_eph = EPH[4]
        # if rok_eph != rok_nav or mes_eph != mes_nav or den_eph != den_nav:
        #    break
        #        print rok_eph, den_eph
        # X_eph = float(EPH[0][1]) * 1000
        # Y_eph = float(EPH[0][2]) * 1000
        # Z_eph = float(EPH[0][3]) * 1000
        # cd_eph = str(EPH[0][0])
        # cas_eph = EPH[1]
        # Rx = X_eph - poloha[0]
        # Ry = Y_eph - poloha[1]
        # Rz = Z_eph - poloha[2]
        # vektor = math.sqrt(math.pow(Rx, 2) + math.pow(Ry, 2) + math.pow(Rz, 2))
        # if vektor > vektor_max:
        #   vektor_max = vektor
        # if vektor < vektor_min:
        #    vektor_min = vektor

        # zapis = cd_eph + "\t" + str(cd_nav) + "\t" + str(cas_eph) + "\t" + str(cas_nav) + "\t" + str(
        #    Dt) + "\t" + str(X_eph) + "\t" + str(Y_eph) + "\t" + str(Z_eph) + "\t"
        # zapis1 = zapis + str(poloha[0]) + "\t" + str(poloha[1]) + "\t" + str(poloha[2]) + "\t"
        # zapis2 = zapis1 + str(Rx) + "\t" + str(Ry) + "\t" + str(Rz) + "\t" + str(vektor) + "\n"
        # xls_porovnanie.write(zapis2)  # zapis vyslednych hodnot xls
        # BLHd = modul2.fXYZ_to_LatLonH(poloha[0], poloha[1], poloha[2])
        zapis_txt = cd_nav + " " + str(cas_nav + Dt) + " " + str(cas_nav) + " " + str(
            Dt) + " " + str(poloha[0]) + " " + str(poloha[1]) + " " + str(
            poloha[2]) + "\n"  # zapis vyslednych hodnot txt
        txt_poloha_d.write(zapis_txt)
        #  vypocet elevacneho uhlu
        Xd = poloha[0]
        Yd = poloha[1]
        Zd = poloha[2]
        dx = Xd - X
        dy = Yd - Y
        dz = Zd - Z
        l1 = math.sqrt(math.pow(dx, 2) + math.pow(dy, 2))  # priemet druzice do vysky ref.stanice
        l2 = math.sqrt(math.pow(dx, 2) + math.pow(dy, 2) + math.pow(dz, 2))
        a_rad = math.acos(l1 / l2)
        a = math.degrees(a_rad)
        #data.append(cd_nav)
        #data.append(cas_nav)
        #data.append(Dt)
        #data.append(Xd)
        #data.append(Yd)
        #data.append(Zd)
        #data.append(a_rad)
        #data_all.append(data)
        if a >= 5 and a <= 30:
            data.append(cd_nav)
            data.append(cas_nav)
            data.append(Dt)
            data.append(Xd)
            data.append(Yd)
            data.append(Zd)
            data.append(a_rad)
            data_all.append(data)
        Dt = Dt + interval
        # print k
    k = k + 1
xls_porovnanie.close()
txt_poloha_d.close()
##print vektor_max
##print vektor_min
##print data_all[0]
# ##########################################################

len_hlavicka_obs = len(hlavicka_obs)
telo = subor_obs[len_hlavicka_obs:]
i = 0
signalS1 = "S1"
signalS2 = "S2"
index_signalS1 = typ_observ.index(signalS1)
index_signalS2 = typ_observ.index(signalS2)
dlzka_zaznam = 16  # dlzka zaznamu je 16 >> na riadku je 5 zaznamov s celkovou dlzkou 80
index_obsS1 = index_signalS1 * dlzka_zaznam
index_obsS2 = index_signalS2 * dlzka_zaznam
axisXelevUhol = []
axisXsinElevUhol = []
axisXcas = []
axisYSNR1 = []
axisYSNR2 = []
# ###########################################################
zoznam_suradnic_xyz = []
zoznam_suradnic_blh = []
zoznam_xyz = [nazov, X, Y, Z]
zoznam_blh = [nazov, L, B, H]
zoznam_suradnic_xyz.append(zoznam_xyz)
zoznam_suradnic_blh.append(zoznam_blh)
r = 6378.135    # polomer zeme v km
lambda1 = 19.0  # vlnova dlzka L1 v cm
lambda2 = 24.4  # vlnova dlzka L2 v cm
while len(telo) > 0:
    data_obs = modul2.fObservacie(telo, pocet_kan)
    for j in range(len(data_all)):
        sur_ort = []
        cas_nav = data_all[j][1]
        Dt = data_all[j][2]
        cas_obs = data_obs[0][3]
        cas = cas_nav + Dt
        if cas != cas_obs:
            continue
        zoznam_druzic = data_obs[0][4]
        cd_nav = data_all[j][0][1:]
        Xd1 = data_all[j][3]
        Yd1 = data_all[j][4]
        Zd1 = data_all[j][5]
        ##print cd_nav, cas_nav
        elev_uhol_rad = data_all[j][-1]
        elev_uhol_stupne = math.degrees(elev_uhol_rad)
        sin_elev_uhol = math.sin(elev_uhol_rad)
        # print sin_elev_uhol
        BLH = modul2.fXYZ_to_LatLonH(Xd1, Yd1, Zd1)
        Bd = BLH[0]
        Ld = BLH[1]
        Hd = BLH[2]
        #Hd_geoid = Hd - modul2.fhel_to_geoid(Bd, Ld)
        Bd_rad = math.radians(Bd)  # Latitude      Zemepisna sirka fi   y
        Ld_rad = math.radians(Ld)  # Longitude     Zemepisna dlzka lambda  x
        B_rad = math.radians(B)
        L_rad = math.radians(L)
        Azimut2 = modul2.fCalculateAzimuth(X, Y, Xd1, Yd1)
        n = 0
        m = 3
        for c in range(len(zoznam_druzic) / 3):
            cd_obs = zoznam_druzic[n:m]
            n = n + 3
            m = m + 3
            if cd_obs != cd_nav:
                continue
                # break
            observacia = data_obs[1][c]
            zaznam_obsS1 = observacia[index_obsS1:index_obsS1 + dlzka_zaznam].split()
            zaznam_obsS2 = observacia[index_obsS2:index_obsS2 + dlzka_zaznam].split()
            # print zaznam_obs
            if len(zaznam_obsS1) > 0:
                for cislo in zaznam_obsS1:
                    # print cislo
                    if len(cislo) < 2:
                        zaznam_obsS1.remove(cislo)
            elif len(zaznam_obsS1) == 0:
                break
            snr1 = float(zaznam_obsS1[0])
            snr_lin1 = math.pow(10, (snr1 / 20))
            # print snr1
            if len(zaznam_obsS2) > 0:
                for cislo in zaznam_obsS1:
                    # print cislo
                    if len(cislo) < 2:
                        zaznam_obsS2.remove(cislo)
            elif len(zaznam_obsS2) == 0:
                break
            snr2 = float(zaznam_obsS2[0])
            axisXcas.append(cas / 60 / 60)
            axisXelevUhol.append(math.degrees(elev_uhol_rad))
            axisXsinElevUhol.append(sin_elev_uhol)
            axisYSNR1.append(snr1)
            axisYSNR2.append(snr2)
            snr_lin2 = math.pow(10, (snr2 / 20))
            h = pilier + float(ant_HEN[0])
            psi1 = ((2*pi)/lambda1)*2*h*math.sin(elev_uhol_rad)
            zapis = cd_obs.replace(".", ",") + "\t" + str(cas_obs).replace(".", ",") + "\t" + \
                    str(sin_elev_uhol).replace(".", ",") + "\t" + str(elev_uhol_stupne).replace(".", ",") + \
                    "\t" + str(snr1).replace(".", ",") + "\t" + str(snr_lin1).replace(".", ",") + \
                    "\t" + str(snr2).replace(".", ",") + "\t" + str(snr_lin2).replace(".", ",") + \
                    "\t" + str(Xd1).replace(".", ",") + "\t" + str(Yd1).replace(".", ",") + "\t" + str(Zd1).replace(".",",") + \
                    "\t" + str(Bd).replace(".", ",") + "\t" + str(Ld).replace(".", ",") + "\t" + str(Hd).replace(".",",") + \
                    "\t" + str(Azimut2).replace(".",",") + "\t" + str(psi1).replace(".",",") + "\n"

            xls_obs.write(zapis)
            zoznam_xyz = [cas, Xd1, Yd1, Zd1]
            zoznam_blh = [cas, Bd, Ld, Hd, cd_obs, Azimut2]
            zoznam_suradnic_xyz.append(zoznam_xyz)
            zoznam_suradnic_blh.append(zoznam_blh)
            break
xls_obs.close()
kruh = kruh.fkruznica(B, L, nazov)
StanicaDruzica = SHP.fLineShp(zoznam_suradnic_blh, nazov)
orezanie = priesecnik.fIntersect(kruh, zoznam_suradnic_blh, nazov)
interval_a = [[80, 100], [160,215]]
az_zoznam = []
for i in range(1,len(orezanie)):
    az = orezanie[i][4]

    if (az > interval_a[0][0] and az < interval_a[0][1]) or (az > interval_a[1][0] and az < interval_a[1][1]):
        az_zoznam.append(orezanie[i])
        print az
az_zoznam.insert(0,[L, B])
linia = SHP.fLineClipShp(az_zoznam, nazov)