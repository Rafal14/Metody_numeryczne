# -*-coding: utf-8 -*-

# Program do wyznaczania wielomianu interpolacyjego
# w postaci Lagrande'a, Newtona i z warunków układu interpolacyjnego

import numpy as np
import matplotlib.pyplot as pt


def intepolacjaLagrande(wekt_x, wekt_y, N, pkt_x):
    '''Oblicza wartość funkcji wielomianu interpolacyjnego w postaci Lagrande'''
    y_k = 0.0

    for k in range(0,N):
        L_k    = 1.0
        for i in range(0,N):
            if i != k:
                L_k = L_k * ( (pkt_x - wekt_x[i])/(wekt_x[k] - wekt_x[i]) )

        y_k = y_k + wekt_y[k] * L_k
    return y_k




def ilorazyRoznicowe(x, y):
    '''Oblicza wartości kolejnych ilorazów różnicowych w met Newtona'''
    # lWez - liczba wezłów
    lWez = len(x)
    B = np.zeros(shape=(lWez, lWez), dtype=float)

    wektor = []

    # wpisanie do pierwszej kolumny macierzy argumentów x
    for licz in range(0,lWez):
        B[licz][0] = y[licz]

    # obliczenie ilorazów różnicowych
    for kol in range(1,lWez):
        for wier in range(kol, lWez):
            B[wier][kol] = (B[wier][kol-1] - B[wier-1][kol-1])/(x[wier] - x[wier-kol])

    for r in range(0,lWez):
        wektor.append(B[r][r])

    return wektor



def interpolacjaNewton( wektorX, ilorazRoz, punkt_x):
    '''Oblicza wartość funkcji wielomianu interpolacyjnego w postaci Newtona'''
    wartWielom = ilorazRoz[0]
    liczWez = len(wektorX) - 1
    roznica = 1
    for licznik in range (1, liczWez):
        wartWielom = wartWielom + ( ( punkt_x - wektorX[licznik-1] ) * roznica ) * ilorazRoz[licznik]
        roznica = ( punkt_x - wektorX[licznik-1] ) * roznica

    return wartWielom



def zamien_wiersze(macierz,w1,w2):
    '''Zamienia wiersze w macierzy A i zwraca zmienioną macierz'''
    l = w1
    while l < (len(macierz)+1):
        pop = macierz[w1][l]
        macierz[w1][l] = macierz[w2][l]
        macierz[w2][l] = pop
        l +=1
    return macierz

def zamien_kolumny(macierz,k1,k2):
    '''Zamienia kolumny w macierzy A i zwraca zmienioną macierz'''
    l = k1
    while l < (len(macierz)):
        pop = macierz[l][k1]
        macierz[l][k1] = macierz[l][k2]
        macierz[l][k2] = pop
        l +=1
    return macierz



def WyznaczWektoraWspol ( wek_X, wek_Y):
    '''Wyznacza wektor współ potrzebny do uzyskania wielom interpret z ukł warunków interpolacyjnych'''
    eq_no= len(wek_X)

    # utwórz macierz Vandermonde'a
    A = np.vander(wektor_x, N=eq_no)

    varray = np.asarray(wek_Y, dtype=float)

    b = np.array([varray])
    # Rozwiązanie

    # Do macierzy A dołącz kolumnę - wektor b
    A1 = np.concatenate((A, b.T), axis=1)

    # znajdź największy współczynnik względem 1. kolumny
    rozmiar = eq_no
    szer = eq_no + 1
    licznik = 0
    ind_wier= 0
    ind_kol = 0

    # wektor zmian kolumn
    historia_kol = []
    # wektor zmian wierszy
    historia_wier = []

    for licznik in range(0, rozmiar - 1):
        pop = A1[licznik][licznik]
        ind_wier = licznik
        ind_kol = licznik

        # znajdź współczynnik o największym module w macierzy A
        for r in range(licznik, rozmiar):
            for c in range(licznik, rozmiar):
                if (pop < abs(A1[r][c])):
                    pop = abs(A1[r][c])
                    ind_wier = r
                    ind_kol = c

        if (ind_wier != licznik):
            zamien_wiersze(A1, licznik, ind_wier)
            historia_wier.append(licznik)
            historia_wier.append(ind_wier)
            if (ind_kol != licznik):
                zamien_kolumny(A1, licznik, ind_kol)
                historia_kol.append(licznik)
                historia_kol.append(ind_kol)
        else:
            if (ind_kol != licznik):
                zamien_kolumny(A1, licznik, ind_kol)
                historia_kol.append(licznik)
                historia_kol.append(ind_kol)

        # eliminacja współczynnika dla x_{i} bez pierwszego wiersza
        for wier in range(licznik + 1, rozmiar):
            krot = A1[wier][licznik] / A1[licznik][licznik]
            for kol in range(0, szer):
                A1[wier][kol] = A1[wier][kol] - (A1[licznik][kol] * krot)

    # wartość ostatniego elementu
    A1[rozmiar - 1][szer - 1] = A1[rozmiar - 1][szer - 1] / A1[rozmiar - 1][szer - 2]
    A1[rozmiar - 1][szer - 2] = 1

    # Podstawienie wsteczne
    cols = eq_no
    rows = eq_no - 2
    count = eq_no - 2

    while cols > 0:
        while rows > -1:
            A1[rows][cols] = A1[rows][cols] - (A1[rows][cols - 1] * A1[rows + 1][cols])
            A1[rows][cols - 1] = 0
            if (rows == 0):
                A1[count][cols] = A1[count][cols] / A1[count][count]
                A1[count][count] = 1
                count -= 1
            rows -= 1
        cols -= 1
        rows = count

    # wyświetl wynik
    # macierz współczynników
    a = []
    for co in range(0, eq_no):
        a.append(A1[eq_no-co-1][eq_no])
    return a



def interpolacjaUkladuWarunkow (wektor_a, punkt):
    '''Oblicza wartość funkcji wielomianu interpolacyjnego wyznaczonego z ukl warunkow interpol'''
    roz = len(wektor_a)

    WartWiel = 0.0
    for p in range (0, roz):
        pot = punkt ** p
        WartWiel = ( wektor_a[p] * pot ) + WartWiel

    return WartWiel







#interpolacja
print "Program do wyznaczania wielomianu interpolacyjengo"
print "trzema metodami\n"

print "Wybierz typ testu:"
print "\n1 - Interpolacja funkcji exp(x) - węzły równoległe"
print "2 - Interpolacja funkcji 1/(1 + x**2) - węzły równoległe"
print "3 - Interpolacja funkcji exp(x) - przesuniete zera wielo. Czebyszewa I rodzaju"
print "4 - Interpolacja funkcji 1/(1 + x**2) - przesuniete zera wielo. Czebyszewa I rodzaju"
test = int(raw_input("\nTwój wybór: "))


print "Podaj wartości a"
wart_a = float(raw_input("a  = "))
print "Podaj wartości b"
wart_b = float(raw_input("b  = "))

# liczba wezlów
wezly = 25

wektor_x = np.zeros(shape=wezly, dtype=float)
wektor_y = np.zeros(shape=wezly, dtype=float)



if test == 1 or test == 2:
    for lp in range(0,wezly):
        wektor_x[lp] = wart_a + (( wart_b - wart_a )/(wezly-1) ) * lp
elif test == 3 or test == 4:
    for lp in range(0,wezly):
        f      = (wart_a - wart_b)/2
        h      = (wart_a+wart_b)/2
        wartPi = np.pi
        wektor_x[lp] = np.cos( ((2.0*lp +2.0)/(2.0*wezly))*wartPi )
        wektor_x[lp] = wektor_x[lp] * f      + h

if test == 1 or test == 3:
    for u in range(0, len(wektor_x)):
        wektor_y[u] = np.exp(wektor_x[u])

elif test == 2 or test == 4:
    for u in range(0, len(wektor_x)):
        wektor_y[u] = 1/( 1 + (wektor_x[u])**2 )


liczba = len(wektor_x) - 1                     # o - mniejsza o 1 niz liczba pomiarow

# wyznaczenie ilorazów różnicowych koniecznych w metodzie Newtona
IlorazyRoznic =ilorazyRoznicowe(wektor_x, wektor_y)


#wyznaczenie wektora wspol a do pkt 3.
wektorA = WyznaczWektoraWspol(wektor_x, wektor_y)



#wykreślenie wykresów wielomianów interpolacyjnych

#zakres osi x
t = np.arange(wart_a, wart_b, 0.1)

#wykres
pt.figure(1)
wyk = pt.subplot(121)
#Wykreślenie pktow
wyk.plot(wektor_x, wektor_y, 'bo', label='pkt')
wyk.grid(True)
wyk.hold(True)
#wielomianu
wyk.plot(t, intepolacjaLagrande(wektor_x, wektor_y, liczba, t), 'r', label='WielomLagrande')
#dodanie legendy
pt.legend()
pt.xlabel('x', fontsize=12, color='black')
pt.ylabel('y')
pt.title('Interpolacja wielomianem w postaci Lagrange', fontsize=12, color='black')



wyk = pt.subplot(122)
#Wykreślenie pktow
wyk.plot(wektor_x, wektor_y, 'ko', label='pkt')
wyk.grid(True)
wyk.hold(True)
wyk.plot(t, interpolacjaNewton(wektor_x, IlorazyRoznic, t), 'k', label='WielomNewtona')

#dodanie legendy
pt.legend()
pt.xlabel('x', fontsize=12, color='black')
pt.ylabel('y')
pt.title('Interpolacja wielomianem w postaci Newtona', fontsize=12, color='black')
pt.show()


#wykres
pt.figure(2)
wyk = pt.subplot(111)
#Wykreślenie pktow
wyk.plot(wektor_x, wektor_y, 'bo', label='pkt')
wyk.grid(True)
wyk.hold(True)
#wielomianu
wyk.plot(t, interpolacjaUkladuWarunkow(wektorA, t), 'r', label='WielomWynZuklWar')
#dodanie legendy
pt.legend()
pt.xlabel('x', fontsize=12, color='black')
pt.ylabel('y')
pt.title('Interpolacja wielomianem', fontsize=12, color='black')
pt.show()