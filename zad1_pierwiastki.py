# -*-coding: utf-8 -*-

# Program do znajdowania pierwiastków równania
# x**3 + x**2 - 3*x - 2 = 0
# w przedziałach [-3, -1] i [-1, 1]
# metodą bisekcji

import numpy as np


def f(u):
    '''Oblicza wartość funkcji wielomianu z zadania'''
    wartosc = u ** 3 + u ** 2 - 3 * u - 2

    return wartosc


print "Program do znajdowania pierwiastków równania\n"
print "x**3  +  x**2  - 3*x  - 2 = 0\n"
print "w przedziałach [-3, -1] i [-1, 1]"
print "metodą bisekcji."

# wartości ograniczeń przedziału
a = -3.0
b = -1.0
c = 1.0


N = int(raw_input("\nPodaj liczbę kroków "))

# obliczenie pierwiastków w przedziale [-3, -1]
# dolna granica przedziału
d = a

# górna granica przedziału
u = b

for licz in range (0,2):
    for k in range(0, N):
        # środek przedziału
        m = (d + u) / 2

        # oblicznie wartości funkcji wielomianu dla m i d
        wart_m = f(m)
        wart_d = f(d)

        # określenie znaku otrzymanych wartości
        znak_wartm = np.sign(wart_m)
        znak_wartd = np.sign(wart_d)

        #sprawdzenie czy wartość funkcji dla argumentu z środka przedziału jest równa 0
        # sprawdzenie czy wartości na końcach przedziałów mają ten sam znak
        # czyli nie funkcja w tym przedziale nie przechodzi przez osi OX
        if znak_wartm == znak_wartd:
            d = m

        # sprawdzenie czy wartości na końcach przedziałów mają przeciwne znaki
        # czyli funkcja w tym przedziale przechodzi przez osi OX
        if znak_wartm == ( znak_wartd * (-1) ):
            u = m

    #przedział, w którym znajduje się pierwiastek
    print "\nWyznaczony przedział, w którym znajduje się pierwiastek"
    print "[d, u] = ", "[", d, ",", u, "]"

    #wyznaczenie pierwiastka
    pierw = (d + u)/2
    wart_p= f(pierw)
    print "\nPierwiastek wynosi ", pierw
    print "f(",pierw,") = ", wart_p

    #zmian granic przedziału dla obliczenia pierwiastka w przedziale [-1,1]
    d = b
    u = c

#wyliczenie miary błędu
bl_1 = (b-a)/(2**N)
bl_2 = (c-b)/(2**N)

print "\nBłąd wyznaczenia pierwiastka w przedziale [-3,-1] wynosi ", bl_1
print "Błąd wyznaczenia pierwiastka w przedziale [-1,1]  wynosi ", bl_2

