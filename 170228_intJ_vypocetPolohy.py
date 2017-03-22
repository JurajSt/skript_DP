# coding=utf-8
import sys
import os
import modul2
import math

## epresne efemeridy
cesta_eph = os.path.join("../data/vstup/cod19212.eph")  # cof19196.eph
if not os.path.exists(cesta_eph):
    print "subor eph neexistuje"
    sys.exit()

## navigacna sprava
cesta_nav = os.path.join("../data/vstup/ganp3060.16n")  # gan60.16n gope0480.14n
if not os.path.exists(cesta_nav):
    print "subor nav neexistuje"
    sys.exit()

vstup_eph = open(cesta_eph, "r")
vstup_nav = open(cesta_nav, "r")
subor_eph = vstup_eph.readlines()
subor_nav = vstup_nav.readlines()
vstup_eph.close()
vstup_nav.close()

pocetEph = 0  # aby bol nazov vystupu totozny z
for znak in cesta_eph:  # vstupnymi subormi
    if "/" in znak:
        pocetEph = pocetEph + 1
nazovEph = cesta_eph.split("/")[pocetEph].split(".")[0]
pocetNav = 0
for znak in cesta_nav:
    if "/" in znak:
        pocetNav = pocetNav + 1
nazovNav = cesta_nav.split("/")[pocetNav].split(".")[0]
## xls vystup
koncovka_xls = ".xls"
nazov_xls = nazovEph + "_" + nazovNav + koncovka_xls
cesta_vystup_xls = "../data/vystup/" + nazov_xls
data_vystup_xls = os.path.join(cesta_vystup_xls)  # vystup xls
xls = open(data_vystup_xls, "w")
xls.write("CD_eph\tCD_nav\tcas_eph\tcas_nav\tDt\tXeph\tYeph\tZeph\tXvyp\tYvyp\tZvyp\tXroz\tYroz\tZroz\tVektor\n")
## txt vystup
koncovka_txt = ".txt"
nazov_txt = nazovNav + "_poloha_druzice"
cesta_vystup_txt = "../data/vystup/" + nazov_txt + koncovka_txt
data_vystup_txt = os.path.join(cesta_vystup_txt)  # vystup txt
txt = open(data_vystup_txt, "w")
txt.write("Cislo_druzice cas_eph cas_nav Dt X_vyp Y_vyp Z_vyp\n")
## txt vystup pre efemeridy
##nazov_eph_txt = nazovNav+"_eph"
##cesta_vystup_eph_txt = "../data/vystup/"+nazov_eph_txt+koncovka_txt
##data_vystup_eph_txt = os.path.join(cesta_vystup_eph_txt)    # vystup txt
##txt_eph = open(data_vystup_eph_txt,"w")
##txt_eph.write("Cislo_druzice cas_eph X_vyp Y_vyp Z_vyp\n")
## ###########################################################
## nacitanie udajov z eph
i = 0
while "*" is not subor_eph[i][0]:
    i = i + 1
telo_eph = subor_eph[i:]
## ##########################################################
## ganp aprox
X = 3929181.900
Y = 1455236.500
Z = 4793653.800
## ##########################################################
## nacitanie udajov z navigacnej spravy
hlavicka = []
for data in subor_nav:
    hlavicka.append(data)
    if "END" in data:
        break
telo_nav = subor_nav[len(hlavicka):]
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
    # kontrola = modul2.fEph_zapis(telo_eph, nazovNav, koncovka_txt)
    while Dt < (platnost_spravy + dodatok):
        data = []
        poloha = modul2.fvypocet_poloha(sprava, Dt)  # vypocet polohy druzice v zadanom case Dt

        EPH = modul2.fData_EPH(telo_eph, cd_nav, (cas_nav + Dt))

        if EPH == "konec zaznamu":
            break
        rok_eph = int(str(EPH[2])[2:])
        mes_eph = EPH[3]
        den_eph = EPH[4]
        if rok_eph != rok_nav or mes_eph != mes_nav or den_eph != den_nav:
            break
            ##        print rok_eph, den_eph
        X_eph = float(EPH[0][1]) * 1000
        Y_eph = float(EPH[0][2]) * 1000
        Z_eph = float(EPH[0][3]) * 1000
        cd_eph = str(EPH[0][0])
        cas_eph = EPH[1]
        Rx = X_eph - poloha[0]
        Ry = Y_eph - poloha[1]
        Rz = Z_eph - poloha[2]
        vektor = math.sqrt(math.pow(Rx, 2) + math.pow(Ry, 2) + math.pow(Rz, 2))
        if vektor > vektor_max:
            vektor_max = vektor
        if vektor < vektor_min:
            vektor_min = vektor

        zapis = cd_eph + "\t" + str(cd_nav) + "\t" + str(cas_eph) + "\t" + str(cas_nav) + "\t" + str(
            Dt) + "\t" + str(X_eph) + "\t" + str(Y_eph) + "\t" + str(Z_eph) + "\t"
        zapis1 = zapis + str(poloha[0]) + "\t" + str(poloha[1]) + "\t" + str(poloha[2]) + "\t"
        zapis2 = zapis1 + str(Rx) + "\t" + str(Ry) + "\t" + str(Rz) + "\t" + str(vektor) + "\n"
        xls.write(zapis2)  # zapis vyslednych hodnot xls
        zapis_txt = cd_eph + " " + str(cas_eph) + " " + str(cas_nav) + " " + str(
            Dt) + " " + str(poloha[0]) + " " + str(poloha[1]) + " " + str(
            poloha[2]) + "\n"  # zapis vyslednych hodnot txt
        txt.write(zapis_txt)

        ##  vypocet elevacneho uhlu
        Xd = poloha[0]
        Yd = poloha[1]
        Zd = poloha[2]
        dx = Xd - X
        dy = Yd - Y
        dz = Zd - Z
        l1 = math.sqrt(math.pow(X, 2) + math.pow(Y, 2) + math.pow(Z, 2))
        l2 = math.sqrt(math.pow(Xd, 2) + math.pow(Yd, 2) + math.pow(Zd, 2))
        cosa = (Xd * (-X) + Yd * (-Y) + Zd * (-Z)) / (l1 * l2)
        a_rad = math.acos(cosa)
        a_stup = math.degrees(a_rad)
        a = a_stup - 90
        data.append(cd_nav)
        data.append(cas_nav)
        data.append(Dt)
        data.append(Xd)
        data.append(Yd)
        data.append(Zd)
        data.append(a)
        data_all.append(data)
        Dt = Dt + 900
    k = k + 1
xls.close()
txt.close()
##txt_eph.close()
##print vektor_max
##print vektor_min
##print data_all[0]

# ##########################################################
# observacne data
cesta_obs = os.path.join("../data/vstup/ganp3060.16o")
if not os.path.exists(cesta_obs):
    print "subor neexistuje"
    sys.exit()
vstup = open(cesta_obs, "r")
subor_obs = vstup.readlines()
vstup.close()

# vystu xls pre observacne data
pocetObs = 0  # aby bol nazov vystupu totozny z
for znak in cesta_obs:  # vstupnymi subormi
    if "/" in znak:
        pocetObs = pocetObs + 1
nazovObs = cesta_eph.split("/")[pocetEph].split(".")[0]
nazov_obs = nazovObs+ "_linElev" + koncovka_xls
cesta_vystup_xls = "../data/vystup/" + nazov_obs
data_vystup_xls = os.path.join(cesta_vystup_xls)  # vystup xls
xls = open(data_vystup_xls, "w")
xls.write("CD\tcas\tSRN\tsin_elevUhol\n")

# nacitanie udajov z hlavicky observacii
hlavicka = []
for data in subor_obs:
    hlavicka.append(data)
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
    for data in hlavicka:
        if "APPROX POSITION XYZ" in data:
            appPosXYZ = data[0:(len(data) - (len(TappPosXYZ)+1))]
        if "ANTENNA: DELTA H/E/N" in data:
            antHEN = data[0:(len(data) - (len(Tant)+1))]
        if "INTERVAL" in data:
            interval = float(data[0:(len(data) - (len(Tint)+1))])
        if "# / TYPES OF OBSERV" in data:
            typObs = data[0:(len(data) - (len(TObs)+1))]
            ts = typObs + " "
            ts1 = ts1 + ts
        if "TIME" in data:
            firstObs = data[0:(len(data) - (len(TfObs)+1))]

first_observ = firstObs.split()
appPos_XYZ = appPosXYZ.split()
ant_HEN = antHEN.split()
typ_observ = ts1.split()
pocet_kan = int(typ_observ[0])
typ_observ.remove(typ_observ[0])
#print "pocet kanalov:", pocet_kan
#print "suradnice XYZ:", appPos_XYZ
#print "pocet a signaly:", typ_observ
#print "prva observacia:", first_observ
#print "interval observacii", interval
#print "vysky anteny", ant_HEN

len_hlavicka = len(hlavicka)
telo = subor_obs[len_hlavicka:]
i = 0
signalS1 = "S1"
index_signal = typ_observ.index(signalS1)
dlzka_zaznam = 16  # dlzka zaznamu je 16 >> na riadku je 5 zaznamov s celkovou dlzkou 80
index_obs = index_signal* dlzka_zaznam
while len(telo)>0:
    data_obs = modul2.fObservacie(telo, pocet_kan)
    for j in range(len(data_all)):
        cd_nav = data_all[j][0][1:]
        cas_nav = data_all[j][1]
        Dt = data_all[j][2]
        cas = cas_nav+Dt
        elev_uhol = data_all[j][-1]
        sin_elev_uhol = math.sin(elev_uhol)
        cas_obs = data_obs[0][3]
        zoznam_druzic = data_obs[0][4]
        if cas != cas_obs:
            continue
        n = 0
        m = 3
        for c in range (len(zoznam_druzic)/3):
            cd_obs = zoznam_druzic[n:m]
            n = n+3
            m = m+3
            if cd_obs != cd_nav:
                continue
                #break
            observacia = data_obs[1][c]
            zaznam_obs = observacia[index_obs:index_obs+dlzka_zaznam].split()
            if len(zaznam_obs)>0:
                for cislo in zaznam_obs:
                    if len(cislo)<2:
                        zaznam_obs.remove(cislo)
            snr1 = float(zaznam_obs[0])
            snr_lin1 = math.pow(10,(snr1/20))
            zapis = cd_obs + "\t" + str(cas_obs) + "\t" + str(snr_lin1) + "\t" + str(sin_elev_uhol) + "\n"
            xls.write(zapis)
            break
xls.close()


