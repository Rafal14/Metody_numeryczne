#-*-coding: utf-8 -*-

# Program do aproksymacji średniokwadratowej
# za pomocą wielomianów
import numpy as np
import matplotlib.pyplot as pt

def Fi(x, i):
    '''Oblicza wartość funkcji x**i'''
    return x**i

def mnozMacierzy (macierzA, macierzB, wym):
    '''Funkcja mnoży macierze kwadratowe'''
    macierz = np.zeros(shape=[wym,wym], dtype=float)
    for i in range(wym):
        for j in range(wym):
            for k in range(wym):
                macierz[i][j] = macierz[i][j] + macierzA[i][k] * macierzB[k][j]
    return macierz


def mnozMacierzyPrzezWekt (macierzA, wektor, wym):
    '''Funkcja mnoży macierz kwadratową przez wektor'''
    wekt = np.zeros(shape=wym, dtype=float)
    for i in range(wym):
        for j in range(wym):
            wekt[i] = wekt[i] + macierzA[i][j] * wektor[j]
    return wekt


def funkcjaExp (wartosc):
    '''Funkcja oblicza wartość funkcji exp()'''
    return np.exp(wartosc)


def funkcja1_x2 (wartosc):
    '''Funkcja oblicza wartość funkcji 1/(x**2 +1)'''
    return 1/( 1 + wartosc**2 )


def zamien_wiersze(macierz,w1,w2):
    '''Zamienia wiersze w macierzy A i zwraca zmienioną macierz'''
    l = w1
    while l < (len(macierz)+1):
        pop = macierz[w1][l]
        macierz[w1][l] = macierz[w2][l]
        macierz[w2][l] = pop
        l +=1
    print " Zamieniono wiersze: ", w1, " ", w2
    return macierz

def zamien_kolumny(macierz,k1,k2):
    '''Zamienia kolumny w macierzy A i zwraca zmienioną macierz'''
    l = k1
    while l < (len(macierz)):
        pop = macierz[l][k1]
        macierz[l][k1] = macierz[l][k2]
        macierz[l][k2] = pop
        l +=1
    print " Zamieniono kolumny ", k1, " ", k2
    return macierz


def funkcjaWielm (wspolczynniki, wartx):
    '''Oblicza wartość funkcji wielomianowej dla zadanego argumentu'''
    wart_wiel = 0.0
    dlugosc = len(wspolczynniki)
    for ind in range(0,dlugosc):
        wart_wiel = wart_wiel + Fi(wartx, ind) * wspolczynniki[ind]
    return wart_wiel


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


def T_k (k, x):
    '''Zwraca odpowiednia wartość T_k'''
    if   k == 0:
        return 1
    elif k == 1:
        return x
    else:
        #utworz wektor wartosci
        r = k+1
        wartosci = np.zeros(shape=r, dtype=float)
        wartosci[0] = 1.0
        wartosci[1] = x
        for i in range(2,r):
            wartosci[i] =  (2 * x * wartosci[i-1]) - wartosci[i-2]
        return wartosci

def WartPrzesPrzesWielCzebysz (k, wartA, wartB, x):
    '''Zwraca wartość dla przesuniętego przeskalowanego wielmianu Czebyszewa I rodzaju'''
    if k < 2:
        wart_zwr = T_k(k,x) * ( ((2*x)/(wartB-wartA)) + ((wartA+wartB)/(wartA-wartB)) )
        return wart_zwr
    else:
        tab = T_k(k,x)
        wart_zwr = tab[k] * (((2 * x) / (wartB - wartA)) + ((wartA + wartB) / (wartA - wartB)))
        return wart_zwr



#wielomian optymalny
def fi_opt(j, wekt):
    '''Zwraca wartość dla wielomianu ortogonalnego'''
    if j < 0:
        return 0
    elif j == 0:
        return 1
    else:
        wynik = 1.0

        r = len(wekt)

        fi_pp = 1
        fi_p  = 1

        alfa_l = 1.0
        alfa_m = 1.0
        beta_l = 1.0
        beta_m = 1.0

        for i in range(1,j+1):
            for p in range(0,r):
                # licznik i mianownik współczynnika alfa
                fi_p = fi_opt(i-1,wekt)
                fi2p = fi_p * fi_p
                alfa_l += (wekt[i] * fi2p)
                alfa_m += fi2p
                # licznik i mianownik współczynnika beta
                beta_l += fi2p

                fi_pp = fi_opt(i-2,wekt)
                fi2pp = fi_pp * fi_pp
                beta_m += fi2pp

            #wartość współczynnika alfa
            alfa = alfa_l/alfa_m
            #wartość współczynnika beta
            beta = beta_l/beta_m
            wynik = ((wekt[i] - alfa) * fi_p) - (beta * fi_pp)

        return wynik



def C_k (wektor, wart_y):
    r = len(wart_y)
    wek_licz = []
    wynik = 0.0
    for i in range(0,r):
        for it in range(1,r):
            wynik += (wart_y[i] * fi_opt(it, wektor))
        wek_licz.append(wynik)
    return  wek_licz


def S_k (wektor):
    r = len(wektor)
    wek_mian = []
    wynik = 0.0
    for i in range(0, r):
        for it in range(1,r):
            fi_op = fi_opt(it, wektor)
            wynik += (fi_op  * fi_op)
        wek_mian.append(wynik)
    return wek_mian




print "Program do aproksymacji średniokwadratowej"
print "za pomocą wielomianów"



print "Podaj wartości krańców przedziału"
print "Podaj wartości a"
wart_a = float(raw_input("a  = "))
print "Podaj wartości b"
wart_b = float(raw_input("b  = "))


# liczba równań w układzie
lrownan = 10
v = []

wekt_wezlow = np.zeros(shape=lrownan, dtype=float)

# utwórz wektor wartości funkcji
wektor_f = np.zeros(shape=lrownan, dtype=float)



print "Wybierz typ testu:"
print "\n1 - węzły równoległe."
print "2 - węzły będące przesuniętym przeskalowanym wielom Czebyszewa I rodzaju"

test = int(raw_input("\nTwój wybór: "))
flaga=1


if   test == 1:
    for lp in range(0,lrownan):
        #oblicz wartosci dla wezla rownoodleglego
        wekt_wezlow[lp] = wart_a + (( wart_b - wart_a )/(lrownan-1) ) * lp
    flaga = 0

elif test == 2:
    for lp in range(0, lrownan):
        #przesunietw wezły 16 wielomianu Czebyszewa
        f = (wart_a - wart_b) / 2
        h = (wart_a + wart_b) / 2
        wartPi = np.pi
        wekt_wezlow[lp] = np.cos(((2.0 * lp + 2.0) / (2.0 * lrownan)) * wartPi)
        wekt_wezlow[lp] = wekt_wezlow[lp] * f + h
    flaga = 0
else:
    flaga = 1

print "\n3 - zadanie 1."
print "4 - zadanie 2."

test = int(raw_input("\nTwój wybór: "))

if   test == 3 and flaga == 0:
    for lp in range(0, lrownan):
        wart = Fi(wekt_wezlow[lp], lp)
        v.append(wart)
    wektFi = np.asarray(v, dtype=float)
elif  test == 4 and flaga == 0:
    for lp in range(0, lrownan):
        wart = WartPrzesPrzesWielCzebysz(lp, wart_a, wart_b, wekt_wezlow[lp])
        v.append(wart)
    wektFi = np.asarray(v, dtype=float)
else:
    flaga=1

print "\n5 - Aproksymacja funkcji exp(x)"
print "6 - Aproksymacja funkcji 1/(1 + x**2)"

test = int(raw_input("\nTwój wybór: "))

if test == 5 and flaga == 0:
    for q in range(0, lrownan):
        wektor_f[q] = funkcjaExp(wekt_wezlow[q])

elif test == 6 and flaga == 0:
    for q in range(0, lrownan):
        wektor_f[q] = funkcja1_x2(wekt_wezlow[q])

else:
    flaga =1


if flaga == 0:
    # utwórz macierz Vandermonde'a
    D1 = np.vander(wektFi, N=lrownan)
    D = np.zeros(shape=[lrownan, lrownan], dtype=float)
    licznik = 0
    for k in range(1, lrownan + 1):
        for j in range(0, lrownan):
            D[j][lrownan - k] = D1[j][licznik]
        licznik += 1

    D_t = D.transpose()

    # A = D_t * D
    A = mnozMacierzy(D_t, D, lrownan)

    b = mnozMacierzyPrzezWekt(D_t, wektor_f, lrownan)




    # Rozwiązanie
    # Rozwiązanie
    m = 1
    k = 1
    t = 1

    # eliminacja zmiennych z kolumn
    while t < (lrownan):
        m = t
        # eliminacja współczynników z jednej kolumny
        while m < (lrownan):
            wsp = A[m][m - k] / A[m - k][m - k]
            for l in range(0, lrownan):
                roz = A[m - k][l] * wsp
                A[m][l] = A[m][l] - roz
            b[m] = b[m] - (b[m - k] * wsp)

            m += 1
            k += 1
        k = 1
        t += 1
    # postępowanie odwrotne
    poprz = 2
    poprzedni = 0

    # wartość ostaniego elementu
    b[lrownan - 1] = b[lrownan - 1] / A[lrownan - 1][lrownan - 1]
    A[lrownan - 1][lrownan - 1] = 1

    licz = lrownan - 2
    kol = lrownan - 1
    wier = lrownan - 2

    while kol > 0:
        while wier > -1:
            b[wier] = b[wier] - (A[wier][kol] * b[kol])
            A[wier][kol] = 0
            if (wier == 0):
                b[licz] = b[licz] / A[licz][licz]
                A[licz][licz] = 1
                licz -= 1
            wier -= 1
        kol -= 1
        wier = licz


    #interpolacja
    # wyznaczenie ilorazów różnicowych koniecznych w metodzie Newtona
    IlorazyRoznic = ilorazyRoznicowe(wekt_wezlow, wektor_f)



    # zakres osi x
    t = np.arange(wart_a, wart_b+0.1, 0.1)

    # wykres
    pt.figure(1)
    # Wykreślenie pktow
    pt.plot(wekt_wezlow, wektor_f, 'ko', label='pkt')
    pt.grid(True)
    pt.hold(True)
    # wielomianu
    pt.plot(t, funkcjaWielm(b, t), 'r', label='WielApro')

    pt.title('Aproksymacja i interpolacja wielomianami', fontsize=12, color='black')

    # wielomianu
    pt.plot(t, intepolacjaLagrande(wekt_wezlow, wektor_f, lrownan-1, t), 'b', label='WielLagrande')

    pt.grid(True)
    pt.hold(True)
    pt.plot(t, interpolacjaNewton(wekt_wezlow, IlorazyRoznic, t), 'g', label='WielNewtona')

    # dodanie legendy
    pt.legend()

    pt.show()

else:
    print "Bład. Nieprawidłowy wybór testu."




#wielomian optymalny
b_k    = []
wsp_fi = []

for no in range(1,lrownan):
    w = fi_opt(no, wekt_wezlow)
    li = w * wektor_f[no]
    mi = w * w

    wsp_fi.append(w)
    b_k.append(li/mi)

