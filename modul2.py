import math
import os, sys
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

def fEph_zapis(telo_eph, nazovNav, koncovka_txt):
    nazov_eph_txt = nazovNav + "_eph"
    cesta_vystup_eph_txt = "../data/vystup/" + nazov_eph_txt + koncovka_txt
    data_vystup_eph_txt = os.path.join(cesta_vystup_eph_txt)  # vystup txt
    txt_eph = open(data_vystup_eph_txt, "w")
    txt_eph.write("Cislo_druzice cas_eph X_vyp Y_vyp Z_vyp\n")
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