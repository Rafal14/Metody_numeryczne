#-*-coding: utf-8 -*-

# Program do rozwiązywania układów równań liniowych
# metodą Gaussa bez wyboru elementu głównego
import numpy as np

print "Program do rozwiązywania układu równań liniowych"
print "metodą Gaussa bez wyboru elementu głównego\n"

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
    b = np.array([b1, b2])
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
    b = np.asarray(u, dtype=float)


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
    b = [1.0, 2.0, 3.0, 4.0]
    eq_no=4
    flaga=0

    print "Macierz A\n", A
    print "\nwektor b\n", b

# Rozwiązanie dla układu z macierzą trójprzekątniową z przekątną, która nie jest dominująca
elif test == 4:
    # macierz trójprzekątniowa z przekątną, która nie jest dominująca
    A = np.array([[1.0, 3.0, 0.0, 0.0], [2.0, 3.0, 4.0, 0.0], [0.0, 4.0, -2.0, -1.0], [0.0, 0.0, 8.0, 6.0]])

    # wektor wyrazów wolnych
    b = [1.0, 2.0, 3.0, 4.0]
    eq_no = 4
    flaga = 0

    print "Macierz A\n", A
    print "\nwektor b\n", b


if (flaga==0):
    #Rozwiązanie
    m = 1
    k = 1
    t = 1

    #eliminacja zmiennych z kolumn
    while t < (eq_no):
        m = t
        #eliminacja współczynników z jednej kolumny
        while m < (eq_no):
            wsp = A[m][m-k] / A[m-k][m-k]
            for l in range(0, eq_no):
                roz = A[m-k][l] * wsp
                A[m][l] = A[m][l] - roz

            b[m] = b[m] - (b[m-k] * wsp)

            m += 1
            k += 1
        k = 1
        t += 1

    #postępowanie odwrotne
    poprz=2
    poprzedni=0

    #wartość ostaniego elementu
    b[eq_no-1] = b[eq_no-1] / A[eq_no-1][eq_no-1]
    A[eq_no - 1][eq_no - 1]=1

    licz=eq_no-2
    kol = eq_no-1
    wier = eq_no-2

    while kol > 0:
        while wier > -1:
            b[wier] = b[wier] - (A[wier][kol] * b[kol])
            A[wier][kol]=0
            if (wier == 0):
                b[licz] = b[licz]/A[licz][licz]
                A[licz][licz]=1
                licz -= 1
            wier -= 1
        kol -= 1
        wier = licz

    #wyświetl wynik
    #macierz współczynników
    print "\nRozwiązanie układu równań\n"
    print "Macierz A\n", A
    print "\nWektor rozwiązań:\n", b
else:
    print "Bład. Nieprawidłowy wybór testu."