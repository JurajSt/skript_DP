import random
def FColorLineKML(c):
    # farba linie Hodnoty R, G, B, T generovane automaticky od 0 do 255, kazda druzica bude mat inu farbu
    R = random.randint(0,255)       # cervena
    G = random.randint(0,255)       # zelena random.randint(0,255)
    B = random.randint(0,255)       # modra random.randint(0,255)
    T = 255                         # transparentnost
    HEX = '%02x%02x%02x' % (B, G, R)   # prepocet farby do hex
    T_HEX = '%02x' % (T)
    farba1 = T_HEX+HEX
    return farba1, c
