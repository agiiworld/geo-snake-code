from proj2_funkcje import FunkcjeGeodezyjne
import numpy as np
import math
import argparse
import sys

parser = argparse.ArgumentParser(
    description="""
Witaj w instrukcji pogramu GEOgik <3. Poniżej znajdziesz odpowiedz na wszystkie twoje pytania (I HOPE SO)
----------------------------------------------------
Wpisuj numer bez kropki :)

Opcje dostępne po uruchomieniu programu z menu:

  1 – Odstęp quasigeoidy (interpolacja dwuliniowa na podstawie pliku WAWA.txt)
  Na podstawie podanej szerokości a potem w długosci (STOPNIE DZIESIĘTNE) z "." jako separatorem dostajemy odstęp do elipsoidy od quasigeoisy.
  Program działa tak, ze wprowadzamy jakies współrzedne (na terenie najlepiej Józefosławia) a program znajduje 4 najbliższe punkty i wykonuje interpolacje.
  ,,Produktem" programu są punkty najbliższe w siatce naszemu wybranemu puntkowi a potem odstęp od quasigeoidy w metrach.
  
  2 – Zamiana fi i lam do układu PL-1992
  Zamienia współrzędne fi i lamda na xy w układzie PL1992. Tutaj Południkiem osiowym
  
  3 – Zamiana fi i lam do układu PL-2000
  Zamienia współrzędne fi i lamda na xy w układzie PL2000. Wymagane jest podanie południka osiowego czyli lambda0.
  
  4 – Zamiana kąta ze stopni dziesiętnych na radiany
  Zwykła kowersja miedzy stopniami i radianami.
  
  5 – Zamiana kąta z radianów na stopnie
  Zwykła konwersja miedzy radianami i stopniami. Tutaj dodstajemy stopnie minuty i sekundy.
  
  6 – Transformacja współrzędnych XYZ na fi, lam i h algorytmem Hirvonena
  Transformacja współrzędnych prostokątnych do współrzędnych eliposidalnych fi i lambda.
  
  7 - Zamienia stopnie minuty i sekundy na stopnie dziesiętne

Dane wejściowe wprowadzane są w formacie dziesiętnym (np. 52.098, 21.032).

Jeżeli chcez miesc stopnie dziesietne bo masz kąt w formacie stp,min,sek użyj opcji 7

Model quasigeoidy WAWA.txt musi znajdować się w tym samym folderze co plik programu. WIĘC GO NIE WYWALAJ ANI NIE ZMIENIAJ NAZWY

Uruchom program bez żadnych parametrów, aby przejść do menu.

Użyj --help, aby zobaczyć tę pomoc.
"""
)
args = parser.parse_args()

def wczytaj_model():
    dane = np.genfromtxt("WAWA.txt", delimiter=" ")
    phi = dane[:, 0]
    lam = dane[:, 1]
    h = dane[:, 2]
    return phi, lam, h

def main():
    geodezja = FunkcjeGeodezyjne()

    while True:
        print("\nSiema geogiku !!! ;) <3 ! Wybierz jedną z dostępnych opcji:")
        print("1. Oblicz odstęp quasigeoidy od elipsoidy")
        print("2. Zamiana fi i lam do układu PL-1992")
        print("3. Zamiana fi i lam do układu PL-2000")
        print("4. Konwersja: stopnie na radiany")
        print("5. Konwersja: radiany na stopnie")
        print("6. Zamiana współrzędnych XYZ na fi, lambda, H")
        print("7. Zamiana stopni, minut i sekund na stopnie dziesiętne")
        print(" ")
        print("Jesli chcesz pomocy napisz --help ;)")
        print(" ")
        print("0. Zakończ program")

        wybor = input("Wpisz numer opcji: ")

        if wybor == "1":
            lat = float(input("Podaj szerokość phi w stopniach dziesiętnych: "))
            lon = float(input("Podaj długość lambda w stopniach dziesiętnych: "))
            phi, lam, h = wczytaj_model()
            try:
                Punkty,R1,R2,P = geodezja.interpolacja_dwuliniowa(lat, lon, phi, lam, h)
                print(f"Punkty najbliższe na których wykonano interpolacje to: {Punkty}")
                print(f"Odstęp quasigeoidy od elipspidy: {P:.3f} m")
            except ValueError as e:
                print("Błąd interpolacji:", e)

        elif wybor == "2":
            fi_deg = float(input("Podaj szerokość phi w stopniach dziesiętnych: "))
            lam_deg = float(input("Podaj długość lambda w stopniach dziesiętnych: "))
            fi_rad = fi_deg * np.pi / 180
            lam_rad = lam_deg * np.pi / 180
            lam0 = 19 * np.pi / 180
            xgk, ygk = geodezja.odwzorowanie_gauss_kruger(fi_rad, lam_rad, lam0)
            x92, y92 = geodezja.PL_1992(xgk, ygk)
            print(f"Współrzędne PL-1992: X = {x92:.3f}, Y = {y92:.3f}")

        elif wybor == "3":
            fi_deg = float(input("Podaj szerokość phi w stopniach dziesiętnych: "))
            lam_deg = float(input("Podaj długość lambda w stopniach dziesiętnych: "))
            strefa = int(input("Podaj numer południka zerowego (15, 18, 21 lub 24): "))
            lam0 = strefa * np.pi / 180
            fi_rad = fi_deg * np.pi / 180
            lam_rad = lam_deg * np.pi / 180
            xgk, ygk = geodezja.odwzorowanie_gauss_kruger(fi_rad, lam_rad, lam0)
            x2000, y2000 = geodezja.PL_2000(xgk, ygk, strefa)
            print(f"Współrzędne PL-2000: X = {x2000:.3f}, Y = {y2000:.3f}")

        elif wybor == "4":
            deg = float(input("Podaj liczbe stopni: "))
            min = float(input("Podaj liczbe minut: "))
            sek = float(input("Podaj podaj liczbe sekund: "))
            rad = geodezja.dms2rad(deg,min,sek)
            print(f"{deg,min,sek} stopni = {rad:.6f} rad")

        elif wybor == "5":
            rad = float(input("Podaj kąt w radianach: "))
            d, m, s = geodezja.rad2dms(rad)
            print(f"{rad:.6f} rad = {int(d)}° {int(m)}' {s:.3f}\"")

        elif wybor == "6":
            x = float(input("Podaj wpółrzedną X "))
            y = float(input("Podaj wpółrzedną Y "))
            z = float(input("Podaj wpółrzedną Z "))
            phi,lam,h,N = geodezja.hirvonen(x,y,z,)
            phi = phi * 180 / np.pi
            lam = lam * 180 / np.pi
            print(f"Twoje wpółrzędne to fi = {phi:.9f} stopni, lam = {lam:.9f} stopni, wysokość = {h:.3f} m. ")

        elif wybor == "7":
            stp = float(input("Podaj liczbe stopni: "))
            min = float(input("Podaj liczbe minut: "))
            sek = float(input("Podaj liczbe sekund: "))
            wynik = stp + (min/60) + (sek/3600)
            print(f"Twoj kąt w formacie stopni dziesiętnych wynosi:{wynik:.8f} ")


        elif wybor == "0":
            print("Zakończyłeś program. Do zobaczenia geogiku!")
            break

        else:
            print("Nieprawidłowa opcja. Spróbuj ponownie.")
            print("Jesli potrzebujesz szczegółwoej pomocy napisz --help.")

if __name__ == "__main__":
    main()
