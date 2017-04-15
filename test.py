import modul2
import math

X = 3929181.8508
Y = 1455236.5099
Xd = 25875123.7835
Yd = -1148521.8129

xf = 3929181.8508
yf = 1455236.5099
xl = 25875123.7835
yl = -1148521.8129
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

print Azimuth
