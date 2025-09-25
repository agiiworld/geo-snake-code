import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import sys
import math as m
pio.renderers.default = 'browser'


dane = np.genfromtxt("WAWA.txt", delimiter=" ")
phi = dane[:,0] #stopnie dziesietne
lam = dane[:,1] #stopnie dziesietne
h = dane[:,2] #metry

#phi0, lam0 = 52.098, 21.032 #do testowania
#h to jest ten odstep wlasnie taka jakby wysokosc
# def interpolacja_dwuliniowa(phi0,lam0,phi,lam,h):
#     Punkty = []
#     for x,y,Q in zip(phi,lam,h):
#         if np.abs(phi-x)<0.01 and np.abs(lam-y)<0.01:
#             Punkty.append([x,y,Q])
#
#     #print(Punkty)
#     # w liscie mamy[x,y,Q]
#     # to jest odstep elipsoidy od quasigeoidy
#     R1 = Punkty[0][2] + (Punkty[2][2]-Punkty[0][2])/(Punkty[2][0]-Punkty[1][0])*(phi0-Punkty[1][0])
#     R2 = Punkty[1][2] + (Punkty[3][2]-Punkty[1][2])/(Punkty[2][0]-Punkty[1][0])*(phi0-Punkty[1][0])
#     P = R1 + (R2-R1)/(Punkty[1][1]-Punkty[0][1])*(lam0-Punkty[0][1])
#     return(Punkty,R1,R2,P)


#Funkcje geodezyjne zapisane jako klasa


class FunkcjeGeodezyjne:
    def __init__(self, model='grs80'):
        model = model.lower()
        if model in ['grs80', 'wgs84']:
            self.a = 6378137.0
            self.e2 = 0.00669438002290
        elif model == 'kras':
            self.a = 6378245.0
            self.e2 = 0.00669342162296
        else:
            print('Nieznany model')
            sys.exit()

    def interpolacja_dwuliniowa(self, phi0, lam0, fi, lam, h):
        Punkty = []
        for x, y, Q in zip(phi, lam, h):
            if np.abs(phi0 - x) < 0.01 and np.abs(lam0 - y) < 0.01:
                Punkty.append([x, y, Q])

        # print(Punkty)
        # w liscie mamy[x,y,Q]

        R1 = Punkty[0][2] + (Punkty[2][2] - Punkty[0][2]) / (Punkty[2][0] - Punkty[1][0]) * (phi0 - Punkty[1][0])
        R2 = Punkty[1][2] + (Punkty[3][2] - Punkty[1][2]) / (Punkty[2][0] - Punkty[1][0]) * (phi0 - Punkty[1][0])
        P = R1 + (R2 - R1) / (Punkty[1][1] - Punkty[0][1]) * (lam0 - Punkty[0][1])
        return (Punkty, R1, R2, P)

    def dms2rad(self,degres):
        rad = degres/180.0*np.pi
        return rad


    def dms2rad(self, d, mi, s):
        degr = d + mi / 60 + s / 3600
        return degr * m.pi / 180

    def rad2dms(self, radiany):
        degr = radiany * 180 / m.pi
        d = np.trunc(degr)
        mins = (degr - d) * 60
        mi = np.trunc(mins)
        s = (mins - mi) * 60
        return (d, mi, s)

    def find_N(self, fi):
        return self.a / (m.sqrt(1 - self.e2 * (m.sin(fi)) ** 2))

    def find_M(self, fi):
        return (self.a * (1 - self.e2)) / (m.sqrt((1 - self.e2 * (m.sin(fi)) ** 2) ** 3))

    def find_R(self, fi):
        N = self.find_N(fi)
        M = self.find_M(fi)
        return m.sqrt(N * M)

    def hirvonen(self, X, Y, Z):
        p = (X ** 2 + Y ** 2) ** 0.5
        fi = m.atan((Z / p) * (1 - self.e2) ** (-1))

        while True:
            N = self.a / (m.sqrt(1 - self.e2 * (m.sin(fi)) ** 2))
            h = p / m.cos(fi) - N
            fi_next = m.atan((Z / p) * (1 - (N * self.e2) / (N + h)) ** (-1))
            if abs(fi_next - fi) < 0.000001 / 3600 * m.pi / 180:
                break
            fi = fi_next

        lam = m.atan2(Y, X)
        N = self.find_N(fi)
        h = p / m.cos(fi) - N
        return fi, lam, h, N

    def hirvonen_odwrotnie(self, h, fi, lam):
        N = self.find_N(fi)
        x = (N + h) * m.cos(fi) * m.cos(lam)
        y = (N + h) * m.cos(fi) * m.sin(lam)
        z = (N * (1 - self.e2) + h) * m.sin(fi)
        return x, y, z


    def odwzorowanie_gauss_kruger(self, fi, lam, lam0):
        b = self.a * m.sqrt(1 - self.e2)
        e2_prime = (self.a ** 2 - b ** 2) / (b ** 2)
        delta_lambda = lam - lam0
        t = m.tan(fi)
        eta2 = e2_prime * m.cos(fi) ** 2
        N = self.find_N(fi)

        A0 = 1 - (self.e2 / 4) - (3 * self.e2 ** 2 / 64) - (5 * self.e2 ** 3 / 256)
        A2 = (3 / 8) * (self.e2 + self.e2 ** 2 / 4 + 15 * self.e2 ** 3 / 128)
        A4 = (15 / 256) * (self.e2 ** 2 + 3 * self.e2 ** 3 / 4)
        A6 = 35 * self.e2 ** 3 / 3072

        sigma = self.a * (A0 * fi - A2 * m.sin(2 * fi) + A4 * m.sin(4 * fi) - A6 * m.sin(6 * fi))

        x_gk = sigma + (delta_lambda ** 2 / 2) * N * m.sin(fi) * m.cos(fi) * (
            1 + (delta_lambda ** 2 / 12) * m.cos(fi) ** 2 * (5 - t ** 2 + 9 * eta2 + 4 * eta2 ** 2) +
            (delta_lambda ** 4 / 360) * m.cos(fi) ** 4 * (
                61 - 58 * t ** 2 + t ** 4 + 270 * eta2 - 330 * eta2 * t ** 2))

        y_gk = delta_lambda * N * m.cos(fi) * (
            1 + (delta_lambda ** 2 / 6) * m.cos(fi) ** 2 * (1 - t ** 2 + eta2) +
            (delta_lambda ** 4 / 120) * m.cos(fi) ** 4 * (5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2))

        return x_gk, y_gk

def przyblizenie_fi(self, xgk):
    """
    Iteracyjne obliczanie szerokości geodezyjnej fi z długości łuku południka xgk.
    """
    A0 = 1 - (self.e2 / 4) - (3 * self.e2 ** 2 / 64) - (5 * self.e2 ** 3 / 256)
    fi = xgk / (self.a * A0)
    A2 = (3 / 8) * (self.e2 + self.e2 ** 2 / 4 + 15 * self.e2 ** 3 / 128)
    A4 = (15 / 256) * (self.e2 ** 2 + 3 * self.e2 ** 3 / 4)
    A6 = 35 * self.e2 ** 3 / 3072

    while True:
        sigma = self.a * (A0 * fi - A2 * m.sin(2 * fi) + A4 * m.sin(4 * fi) - A6 * m.sin(6 * fi))
        fi_prev = fi
        fi = fi_prev + (xgk - sigma) / (self.a * A0)
        if abs(fi - fi_prev) < (0.000001 / 3600 * m.pi / 180):
            break
    return fi

def przelicz_odwrotnie_gauss_kruger(self, xgk, ygk, lambda0_rad):
    """
    Przelicza współrzędne prostokątne (xgk, ygk) na współrzędne geodezyjne (fi, lam).
    """
    phi1 = self.przyblizenie_fi(xgk)

    N1 = self.a / np.sqrt(1 - self.e2 * np.sin(phi1) ** 2)
    M1 = self.a * (1 - self.e2) / (1 - self.e2 * np.sin(phi1) ** 2) ** 1.5
    t1 = np.tan(phi1)
    eta1 = self.e2 * np.cos(phi1) ** 2 / (1 - self.e2)

    phi = phi1 - (ygk ** 2 * t1) / (2 * M1 * N1) * (
        1 - ygk ** 2 / (12 * N1 ** 2) * (5 + 3 * t1 ** 2 + eta1 - 9 * eta1 * t1 ** 2)
        + ygk ** 4 / (360 * N1 ** 4) * (61 + 90 * t1 ** 2 + 45 * t1 ** 4)
    )

    lambda_value = lambda0_rad + ygk / (N1 * np.cos(phi1)) * (
        1 - ygk ** 2 / (6 * N1 ** 2) * (1 + 2 * t1 ** 2 + eta1)
        + ygk ** 4 / (120 * N1 ** 4) * (5 + 28 * t1 ** 2 + 24 * t1 ** 4 + 6 * eta1 + 8 * eta1 * t1 ** 2)
    )

    return phi, lambda_value

def PL_1992(self, x_gk, y_gk):
    """
    Funkcja cechuje współrzędne płaskie do układu PL-1992.
    """
    m_92 = 0.9993
    x_92 = m_92 * x_gk - 5300000
    y_92 = m_92 * y_gk + 500000
    return x_92, y_92

def PL_2000(self, x_gk, y_gk, nr):
    """
    Funkcja cechuje współrzędne płaskie do układu PL-2000.
    nr – numer strefy (5, 6, 7, 8)
    """
    m_2000 = 0.999923
    x_2000 = m_2000 * x_gk
    y_2000 = m_2000 * y_gk + (nr * 1000000) + 500000
    return x_2000, y_2000

def odcechowanie_z_2000(self, x, y):
    """
    Funkcja usuwa cechowanie z układu PL-2000 i wyznacza lambda0 dla strefy.
    """
    m = 0.999923
    xgk = x / m
    nr = round(y / 1000000)
    if nr == 5:
        lam0 = 15 / 180 * np.pi
    elif nr == 6:
        lam0 = 18 / 180 * np.pi
    elif nr == 7:
        lam0 = 21 / 180 * np.pi
    elif nr == 8:
        lam0 = 24 / 180 * np.pi
    else:
        lam0 = None  # można też rzucić wyjątek
    ygk = (y - 500000 - nr * 1000000) / m
    return xgk, ygk, lam0

def odcechowanie_z_1992(self, X1_1992, Y1_1992):
    """
    Funkcja usuwa cechowanie z układu PL-1992.
    """
    x1gk = (X1_1992 + 5300000) / 0.9993
    y1gk = (Y1_1992 - 500000) / 0.9993
    return x1gk, y1gk

