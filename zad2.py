#-*-coding: utf-8 -*-

# Program do wyznaczania wyznacznika i macierzy odwrotnej do zadanej macierzy
# z wykorzystaniem macierzy LU otrzymanej w wyniku dekompozycji

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

import numpy as np
import scipy.linalg as sp

print "Program do wyznaczania macierzy odwrotnej"
print "z wykorzystaniem macierzy LU otrzymanej\n"
print " w wyniku dekompozycji\n"

print "Wybierz typ testu:"
print "\n1 - Wyznaczenie wyznacznika i macierzy odwrotnej do macierzy A"
print "2 - Wyznaczenie wyznacznika i macierzy odwrotnej do macierzy Vandermonde’a (większy niż 15x15)"
print "3 - Wyznaczenie wyznacznika i macierzy odwrotnej do macierzy trójprzekątniowej z dominującą przekątną"
print "4 - Wyznaczenie wyznacznika i macierzy odwrotnej do macierzy trójprzekątniowej z przekątną, która nie jest dominująca"
test = int(raw_input("\nTwój wybór: "))
flaga=1


if test == 1:
    print "Podaj wartości współczynników macierzy A"
    a11 = float(raw_input("a11 = "))
    a12 = float(raw_input("a12 = "))
    a21 = float(raw_input("a21 = "))
    a22 = float(raw_input("a22 = "))

    #liczba równań
    eq_no = 2

    #macierz współczynników
    A = np.array([[a11, a12], [a21, a22]])
    print "Macierz A\n", A

    flaga=0

#Rozwiązanie dla układu z macierzą Vandermonde'a
elif test == 2:
    # liczba równań w układzie
    eq_no = int(raw_input("Podaj liczbę równań: "))

    #utwórz wektor
    v=[]
    print "\nPodaj wartości wektora do utworzenia macierzy Van der Monde'a"
    for i in range(0,eq_no):
        wart = float(raw_input("Wprowadź wartość: "))
        v.append(wart)
    x = np.asarray(v, dtype=float)

    # utwórz macierz Vandermonde'a
    A = np.vander(x,N=eq_no)
    print "\nMacierz Vandermonde'a:\n", A
    flaga=0

# Rozwiązanie dla układu z macierzą trójprzekątniową z dominującą przekątną
elif test == 3:
    #macierz trójprzekątniowa z dominującą przekątną
    A = np.array([[-5.0, 3.0, 0.0, 0.0], [2.0, 9.0, 3.0, 0.0], [0.0, 4.0, 10.0, -2.0], [0.0, 0.0, 2.0, 7.0]])

    eq_no=4
    flaga=0

    print "Macierz A\n", A

# Rozwiązanie dla układu z macierzą trójprzekątniową z przekątną, która nie jest dominująca
elif test == 4:
    # macierz trójprzekątniowa z przekątną, która nie jest dominująca
    A = np.array([[1.0, 3.0, 0.0, 0.0], [2.0, 3.0, 4.0, 0.0], [0.0, 4.0, -2.0, -1.0], [0.0, 0.0, 8.0, 6.0]])

    eq_no = 4
    flaga = 0

    print "Macierz A\n", A


if (flaga==0):
    #Rozwiązanie

    #dekompozycja LU z wykorzystaniem funkcji, P - macierz jednostkowa, permutacji
    P, L, U = sp.lu(A)

    print "\nmacierz P:"
    print P

    print "\nmacierz L:"
    print L

    print "\nmacierz U:"
    print U

    #wymiary macierzy permutacji
    dlugosc = len(P)
    b=[]
    for dl in range(0,dlugosc):
        b.append(0)

    #odwrotna
    odwA = P

    #oblicz wektory kolumnowe dla macierzy odwrotnej
    for count in range(0, dlugosc):

        # nowa macierz współczynników B
        B = np.dot(L, U)

        #utworz wektor b
        for ind in range(0,dlugosc):
            b[ind] = P[ind][count]

        # Rozwiązanie
        m = 1
        k = 1
        t = 1

        # eliminacja zmiennych z kolumn
        while t < (dlugosc):
            m = t
            # eliminacja współczynników z jednej kolumny
            while m < (eq_no):
                wsp = B[m][m - k] / B[m - k][m - k]
                for l in range(0, eq_no):
                    roz = B[m - k][l] * wsp
                    B[m][l] = B[m][l] - roz

                b[m] = b[m] - (b[m - k] * wsp)

                m += 1
                k += 1
            k = 1
            t += 1

        # postępowanie odwrotne
        poprz = 2
        poprzedni = 0

        # wartość ostaniego elementu
        b[eq_no - 1] = b[eq_no - 1] / B[eq_no - 1][eq_no - 1]
        B[eq_no - 1][eq_no - 1] = 1

        licz = eq_no - 2
        kol = eq_no - 1
        wier = eq_no - 2

        while kol > 0:
            while wier > -1:
                b[wier] = b[wier] - (B[wier][kol] * b[kol])
                B[wier][kol] = 0
                if (wier == 0):
                    b[licz] = b[licz] / B[licz][licz]
                    B[licz][licz] = 1
                    licz -= 1
                wier -= 1
            kol -= 1
            wier = licz

        for no in range(0, dlugosc):
            odwA[no][count] = b[no]

    mod =dlugosc % 2
    licznik=0
    ost= dlugosc-1
    srodek = (dlugosc+1)/2
    if ( mod != 0):
        while licznik!=srodek:
            zamien_kolumny(odwA,licznik,ost-licznik)
            licznik += 1
    else:
        while licznik != (dlugosc/2):
            zamien_kolumny(odwA, licznik, ost - licznik)
            licznik += 1

        # oblicz wyznacznik macierz A

        wyznacznik=U[0][0]
        for r in range(1, dlugosc):
            wyznacznik=wyznacznik * U[r][r]

        # wyświetl wynik
        # macierz współczynników
        print "\nRozwiązanie\n"
        print "Macierz odwrotna od A\n", odwA

        print "\nWyznacznik macierzy A wynosi: ", wyznacznik
else:
    print "Bład. Nieprawidłowy wybór testu."
