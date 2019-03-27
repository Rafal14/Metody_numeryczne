#-*-coding: utf-8 -*-

# Program do rozwiązywania układów równań liniowych
# metodą Gaussa z wyborem elementu głównego

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

print "Program do rozwiązywania układu równań liniowych"
print "metodą Gaussa z wyborem elementu głównego\n"

print "Postać układu równań: Ax = b"

print "Wybierz typ testu:"
print "\n1 - Układ dwóch równań liniowych"
print "2 - Układ z macierzą Vandermonde’a (większy niż 15x15)"
print "3 - Układ z macierzą trójprzekątniową z dominującą przekątną"
print "4 - Układ z macierzą trójprzekątniową z przekątną, która nie jest dominująca"
test = int(raw_input("\nTwój wybór: "))
flaga=1

#Rozwiązanie dla układu dwóch równań liniowych
if test == 1:
    print "Podaj wartości współczynników macierzy A"
    a11 = float(raw_input("a11 = "))
    a12 = float(raw_input("a12 = "))
    a21 = float(raw_input("a21 = "))
    a22 = float(raw_input("a22 = "))
    print "Podaj wartości wyrazów wolnych"
    b1  = float(raw_input("b1  = "))
    b2  = float(raw_input("b2  = "))

    #liczba równań
    eq_no = 2

    #macierz współczynników
    A = np.array([[a11, a12], [a21, a22]])

    #wektor wyrazów wolnych
    b = np.array([[b1, b2]])
    print "Macierz A\n", A
    print "\nWektor wyrazów wolnych\n", b
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

    #utwórz wektor b
    u=[]
    print "\nPodaj wartości wyrazów wolnych"
    for q in range(0, eq_no):
        wol = float(raw_input("Wprowadź wartość: "))
        u.append(wol)
    b = np.array([u])

    # utwórz macierz Vandermonde'a
    A = np.vander(x,N=eq_no)
    print "\nMacierz Vandermonde'a:\n", A
    print "Wektor wyrazów wolnych:\n", b
    flaga=0

# Rozwiązanie dla układu z macierzą trójprzekątniową z dominującą przekątną
elif test == 3:
    #macierz trójprzekątniowa z dominującą przekątną
    A = np.array([[-5.0, 3.0, 0.0, 0.0], [2.0, 9.0, 3.0, 0.0], [0.0, 4.0, 10.0, -2.0], [0.0, 0.0, 2.0, 7.0]])

    #wektor wyrazów wolnych
    b = np.array([[1.0, 2.0, 3.0, 4.0]])
    eq_no=4
    flaga=0

    print "Macierz A\n", A
    print "\nwektor b\n", b

# Rozwiązanie dla układu z macierzą trójprzekątniową z przekątną, która nie jest dominująca
elif test == 4:
    # macierz trójprzekątniowa z przekątną, która nie jest dominująca
    A = np.array([[1.0, 3.0, 0.0, 0.0], [2.0, 3.0, 4.0, 0.0], [0.0, 4.0, -2.0, -1.0], [0.0, 0.0, 8.0, 6.0]])

    # wektor wyrazów wolnych
    b = np.array([[1.0, 2.0, 3.0, 4.0]])
    eq_no = 4
    flaga = 0

    print "Macierz A\n", A
    print "\nwektor b\n", b


if (flaga==0):
    #Rozwiązanie

    #Do macierzy A dołącz kolumnę - wektor b
    A1=np.concatenate((A, b.T), axis=1)

    #znajdź największy współczynnik względem 1. kolumny
    rozmiar=eq_no
    szer=eq_no+1
    licznik = 0
    ind_wier = 0
    ind_kol = 0

    #wektor zmian kolumn
    historia_kol=[]
    # wektor zmian wierszy
    historia_wier=[]

    for licznik in range(0,rozmiar-1):
        pop = A1[licznik][licznik]
        ind_wier = licznik
        ind_kol = licznik

        #znajdź współczynnik o największym module w macierzy A
        for r in range(licznik,rozmiar):
            for c in range(licznik,rozmiar):
                if (pop < abs(A1[r][c])):
                    pop = abs(A1[r][c])
                    ind_wier = r
                    ind_kol  = c

        if (ind_wier!=licznik):
            zamien_wiersze(A1,licznik,ind_wier)
            historia_wier.append(licznik)
            historia_wier.append(ind_wier)
            if (ind_kol!=licznik):
                zamien_kolumny(A1,licznik,ind_kol)
                historia_kol.append(licznik)
                historia_kol.append(ind_kol)
        else:
            if (ind_kol!=licznik):
                zamien_kolumny(A1, licznik, ind_kol)
                historia_kol.append(licznik)
                historia_kol.append(ind_kol)

        #eliminacja współczynnika dla x_{i} bez pierwszego wiersza
        for wier in range(licznik+1, rozmiar):
            krot = A1[wier][licznik] / A1[licznik][licznik]
            for kol in range(0,szer):
                A1[wier][kol] = A1[wier][kol] - (A1[licznik][kol] * krot)

    #wartość ostatniego elementu
    A1[rozmiar-1][szer-1] = A1[rozmiar-1][szer-1] / A1[rozmiar-1][szer-2]
    A1[rozmiar-1][szer-2] = 1

    #Podstawienie wsteczne
    cols=eq_no
    rows=eq_no-2
    count=eq_no-2

    while cols>0:
        while rows>-1:
            A1[rows][cols]=A1[rows][cols]-( A1[rows][cols-1] * A1[rows+1][cols] )
            A1[rows][cols-1]=0
            if (rows==0):
                A1[count][cols]=A1[count][cols] / A1[count][count]
                A1[count][count]=1
                count-=1
            rows -= 1
        cols -= 1
        rows = count

    #wyświetl wynik
    #macierz współczynników
    print "\nRozwiązanie układu równań\n"
    print "Macierz A1\n", A1

else:
    print "Bład. Nieprawidłowy wybór testu."