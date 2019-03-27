#-*-coding: utf-8 -*-

# Program do rozwiązywania układów równań liniowych
# metodą Gaussa z częściowym wyborem elementu głównego


def zamien_wiersze(macierz,w1,w2):
    '''Zamienia wiersze w macierzy A i zwraca zmienioną macierz'''
    l = 0
    while l < (len(macierz)+1):
        pop = macierz[w1][l]
        macierz[w1][l] = macierz[w2][l]
        macierz[w2][l] = pop
        l +=1
    return macierz

import numpy as np

print "Program do rozwiązywania układu równań liniowych"
print "metodą Gaussa z częściowym wyborem elementu głównego\n"

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
    ind = 0

    for j in range(1,rozmiar):
        pop = A1[j-1][j-1]
        for i in range(j,eq_no):
            if ( pop < abs(A1[i][j-1]) ):
                pop = abs(A1[i][j-1])
                ind = i
        #wykonaj zamianę wierszy, tak aby wiersz z największym współczynnikiem był w 1. wierszu
        if ( ind > 0 ):
            zamien_wiersze(A1,j-1,ind)

        #wyznacz krotność
        krot = A1[j][j-1] / A1[j-1][j-1]
        for wier in range(0,szer):
            A1[j][wier] = A1[j][wier] - (A1[j-1][wier] * krot)

    #wyznacz wartość ostatniego elementu macierzy A1
    A1[eq_no-1][eq_no] = A1[eq_no-1][eq_no] / A1[eq_no-1][eq_no-1]
    A1[eq_no-1][szer-2] = 1

    #podstawienie wsteczne
    licz = eq_no - 2
    kol = eq_no - 1
    wier = eq_no - 2
    while kol > 0:
        while wier > -1:
            A1[wier][eq_no] = A1[wier][eq_no] - (A1[wier][kol] * A1[kol][eq_no])
            A1[wier][kol] = 0
            if (wier == 0):
                if (A1[licz][licz] != 0):
                    A1[licz][eq_no] = A1[licz][eq_no] / A1[licz][licz]
                else:
                    print "\nUwaga! Wystąpił problem dzielenia przez 0"
                A1[licz][licz] = 1
                licz -= 1
            wier -= 1
        kol -= 1
        wier = licz

    #wyświetl wynik
    #macierz współczynników
    print "\nRozwiązanie układu równań\n"
    print "Macierz A1\n", A1

else:
    print "Bład. Nieprawidłowy wybór testu."