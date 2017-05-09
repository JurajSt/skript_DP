import sys, math, SHP, modul2

def fIntersect(zoznam_kruh, zoznam_suradnic, n):
    zoznamS = []
    nazov = "priesecnik" + n
    cesta_vystup_txt_priesecnik = "data/vystup/kruh/" + nazov + ".txt"
    txt_priesecnik_blh = open(cesta_vystup_txt_priesecnik, "w")
    txt_priesecnik_blh.write("cas L B\n")
    Ax = zoznam_suradnic[0][1]
    Ay = zoznam_suradnic[0][2]
    A = [Ax, Ay]
    for e in range(len(zoznam_suradnic)-1):
        Bx = zoznam_suradnic[e+1][1]
        By = zoznam_suradnic[e+1][2]
        cas = zoznam_suradnic[e+1][0]
        cd = zoznam_suradnic[e+1][4]
        az = zoznam_suradnic[e+1][5]
        B = [Bx, By]
        for j in range(len(zoznam_kruh)-1):
            Cx = zoznam_kruh[j][0]
            Cy = zoznam_kruh[j][1]
            Dx = zoznam_kruh[j + 1][0]
            Dy = zoznam_kruh[j + 1][1]
            C = [Cx, Cy]
            D = [Dx, Dy]
            S = [0, 0]
            if A[1] == B[1]:  # rovonbeznost z osou x a = 0 pre y a = 1 ale nedostnem dobry priesecnik, riesim nizsie
                a1 = 0
            elif A[0] == B[0]:
                a1 = 1
            else:
                a1 = (B[1] - A[1]) / (B[0] - A[0])
            if C[1] == D[1]:
                a2 = 0
            elif C[0] == D[0]:
                a2 = 1
            else:
                a2 = (D[1] - C[1]) / (D[0] - C[0])
            if a1 == a2:
                print "prusecnik neexistuje"
                sys.exit()
            b1 = A[1] - a1 * A[0]
            b2 = C[1] - a2 * C[0]
            S[0] = (b1 - b2) / (a2 - a1)
            S[1] = a2 * S[0] + b2

            if A[0] == B[0]:  # rovnobeznost z osou y
                S[0] = A[0]
                S[1] = a2 * S[1] + b2
            if C[0] == D[0]:
                S[0] = C[0]
                S[1] = a1 * S[0] + b1
            if (A[0] - S[0]) * (S[0] - B[0]) >= 0:
                a = 1
            # print "lezi mezi body AB"
            else:
                continue
                print "prusecnik nelezi mezi body AB"
            if (C[0] - S[0]) * (S[0] - D[0]) >= 0:
                # print "lezi mezi body CD"
                a = 1
            else:
                continue
                print "prusecnik nelezi mezi body CD"
            azimut = modul2.fCalculateAzimuth(Ax, Ay, S[0], S[1])
            S.append(cas)
            S.append(cd)
            S.append(azimut)

            zapis1 = str(S[2]).replace(".", ",") + " " + str(S[0]).replace(".", ",") + " " + str(S[1]).replace(".", ",") + \
                     str(S[3]).replace(".", ",") + str(S[4]).replace(".", ",") + "\n"
            txt_priesecnik_blh.write(zapis1)
            zoznamS.append(S)
    zoznamS.insert(0, A)
    liniaOrez = SHP.fLineClipShp(zoznamS, nazov)
    txt_priesecnik_blh.close()
    return zoznamS