#-*-coding: utf-8 -*-


import numpy as np
from scipy.misc import derivative
import matplotlib.pyplot as pt

#Definicje funkcji
def f(x):
    '''Zwraca wartość dla funkcji eksponencjalnej'''
    return np.exp(x)

def g(x):
    '''Zwraca wartość dla funkcji exp(-(x**2))'''
    return np.exp((x**2)*(-1))

def j(x):
    '''Zwraca wartość dla funkcji (x**2)*log(x)'''
    logarytm=np.log10(x)
    return (x**2)*logarytm

def k(x):
    '''Zwraca wartość dla funkcji 1/((1 + x**2)**(1/2))'''
    pierwiastek=np.sqrt(1 + (x**2))
    return 1/pierwiastek


print "Przybliżanie pochodnej funkcji e^(x)"
print "Podaj przedział do obliczenia pochodnej"

#określenie przedziału
dolna_gr = float(raw_input("Dolna granica przedziału: "))
gorna_gr = float(raw_input("Górna granica przedziału: "))

#określa otoczenie punktu x0
h = (gorna_gr - dolna_gr)/2

#wyznaczenie środka przedziału
x0 = (dolna_gr + gorna_gr)/2
print "Środek przedziału x0 = ", x0
print "h = ", h

#wartość funkcji f w punkcie x0
wart_f_x0 = f(x0)

z1=dolna_gr-1
z2=gorna_gr+1
#zakres osi x
to = np.arange(z1,z2, 0.1)
t = np.arange(dolna_gr, gorna_gr, 0.1)

#obliczenie pochodnej funkcji f w pkt x0
pochodnaf_Oh = (f(x0+h)-f(x0))/h
pochodnaf_Oh2= (f(x0+h)-f(x0-h))/(h*2)

#wypisanie przybliżenia wartości pochodnej w pkt x0
print "f'(x0) = ", pochodnaf_Oh, "z dok. O(h)"
print "f'(x0) = ", pochodnaf_Oh2, "z dok. O(h**2)"

#wykreślenie stycznych określonych przez pochodną funkcji
# f z dokładnością O(h) i O(h**2)
pt.figure(1)
wyk = pt.subplot(111)
#Wykreślenie charakterystyki funkcji exp(x)
wyk.plot(to, f(to), 'r', label='exp(x)')
wyk.grid(True)
wyk.hold(True)
#Wykreślenie stycznej do funkcji exp(x) z dokładnością O(h)
wyk.plot(t, pochodnaf_Oh*t + wart_f_x0, 'b', label='Przyblizenie pochodnej \nfunkcji exp(x) z dok. O(h)')

#Wykreślenie stycznej do funkcji exp(x) z dokładnością O(h^2)
wyk.plot(t, pochodnaf_Oh2 * t + wart_f_x0, 'k', label='Przyblizenie pochodnej \nfunkcji exp(x) z dok. O(h^2)')

#dodanie legendy
wyk.legend()
pt.xlabel('Czas [s]', fontsize=12, color='black')
pt.ylabel('Wartosc funkcji')
pt.title('Charakterystyka exp(x) i jej pochodne', fontsize=12, color='black')
pt.show()

#===============================================================================
print "\nPrzybliżanie pochodnej funkcji e^(-x**2)"
print "Podaj przedział do obliczenia pochodnej"

#określenie przedziału
dolna_gr = float(raw_input("Dolna granica przedziału: "))
gorna_gr = float(raw_input("Górna granica przedziału: "))

#określa otoczenie punktu x0
h = (gorna_gr - dolna_gr)/2

#wyznaczenie środka przedziału
x0 = (dolna_gr + gorna_gr)/2
print "Środek przedziału x0 = ", x0
print "h = ", h

#wartość funkcji g w punkcie x0
wart_g_x0 = g(x0)

z1=dolna_gr-1
z2=gorna_gr+1
#zakres osi x
to = np.arange(z1,z2, 0.1)
t = np.arange(dolna_gr, gorna_gr, 0.1)

#obliczenie pochodnej funkcji g w pkt x0
pochodnag_Oh = (g(x0+h)-g(x0))/h
pochodnag_Oh2= (g(x0+h)-g(x0-h))/(h*2)

#wypisanie przybliżenia wartości pochodnej w pkt x0
print "g'(x0) = ", pochodnag_Oh, "z dok. O(h)"
print "g'(x0) = ", pochodnag_Oh2, "z dok. O(h**2)"

#wykreślenie stycznych określonych przez pochodną funkcji
# g z dokładnością O(h) i O(h**2)
pt.figure(2)
wyk = pt.subplot(111)
#Wykreślenie charakterystyki funkcji exp(-x**2)
wyk.plot(to, g(to), 'r', label='exp(-x**2)')
wyk.grid(True)
wyk.hold(True)
#Wykreślenie stycznej do funkcji exp(-x**2) z dokładnością O(h)
wyk.plot(t, pochodnag_Oh*t + wart_g_x0, 'b', label='Przyblizenie pochodnej \nfunkcji exp(-x**2) z dok. O(h)')

#Wykreślenie stycznej do funkcji exp(x) z dokładnością O(h^2)
wyk.plot(t, pochodnag_Oh2 * t + wart_g_x0, 'k', label='Przyblizenie pochodnej \nfunkcji exp(-x**2) z dok. O(h^2)')

#dodanie legendy
wyk.legend()
pt.xlabel('Czas [s]', fontsize=12, color='black')
pt.ylabel('Wartosc funkcji')
pt.title('Charakterystyka exp(-x**2) i jej pochodne', fontsize=12, color='black')
pt.show()


#====================================================================
print "\nPrzybliżanie pochodnej funkcji (x^2)log(x)"
print "Podaj przedział do obliczenia pochodnej"

#określenie przedziału
dolna_gr = float(raw_input("Dolna granica przedziału: "))
gorna_gr = float(raw_input("Górna granica przedziału: "))

#określa otoczenie punktu x0
h = (gorna_gr - dolna_gr)/2

#wyznaczenie środka przedziału
x0 = (dolna_gr + gorna_gr)/2
print "Środek przedziału x0 = ", x0
print "h = ", h

#wartość funkcji j w punkcie x0
wart_j_x0 = j(x0)

z1=0.001
z2=gorna_gr+1
#zakres osi x
to = np.arange(z1,z2, 0.1)
t = np.arange(dolna_gr, gorna_gr, 0.1)

#obliczenie pochodnej funkcji j w pkt x0
pochodnaj_Oh = (j(x0+h)-j(x0))/h
pochodnaj_Oh2= (j(x0+h)-j(x0-h))/(h*2)

#wypisanie przybliżenia wartości pochodnej w pkt x0
print "j'(x0) = ", pochodnaj_Oh, "z dok. O(h)"
print "j'(x0) = ", pochodnaj_Oh2, "z dok. O(h**2)"

#wykreślenie stycznych określonych przez pochodną funkcji
# j z dokładnością O(h) i O(h**2)
pt.figure(3)
wyk = pt.subplot(111)
#Wykreślenie charakterystyki funkcji (x^2)log(x)
wyk.plot(to, j(to), 'r', label='(x**2)log(x)')
wyk.grid(True)
wyk.hold(True)
#Wykreślenie stycznej do funkcji (x**2)log(x) z dokładnością O(h)
wyk.plot(t, pochodnaj_Oh*t + wart_j_x0, 'b', label='Przyblizenie pochodnej \nfunkcji (x**2)log(x) z dok. O(h)')

#Wykreślenie stycznej do funkcji (x**2)log(x) z dokładnością O(h^2)
wyk.plot(t, pochodnaj_Oh2 * t + wart_j_x0, 'k', label='Przyblizenie pochodnej \nfunkcji (x**2)log(x) z dok. O(h^2)')

#dodanie legendy
wyk.legend()
pt.xlabel('Czas [s]', fontsize=12, color='black')
pt.ylabel('Wartosc funkcji')
pt.title('Charakterystyka (x**2)log(x) i jej pochodne', fontsize=12, color='black')
pt.show()

#====================================================================
print "\nPrzybliżanie pochodnej funkcji 1/(sqrt(1 + x^2))"
print "Podaj przedział do obliczenia pochodnej"

#określenie przedziału
dolna_gr = float(raw_input("Dolna granica przedziału: "))
gorna_gr = float(raw_input("Górna granica przedziału: "))

#określa otoczenie punktu x0
h = (gorna_gr - dolna_gr)/2

#wyznaczenie środka przedziału
x0 = (dolna_gr + gorna_gr)/2
print "Środek przedziału x0 = ", x0
print "h = ", h

#wartość funkcji k w punkcie x0
wart_k_x0 = k(x0)

z1=dolna_gr-1
z2=gorna_gr+1
#zakres osi x
to = np.arange(z1,z2, 0.1)
t = np.arange(dolna_gr, gorna_gr, 0.1)

#obliczenie pochodnej funkcji k w pkt x0
pochodnak_Oh = (k(x0+h)-k(x0))/h
pochodnak_Oh2= (k(x0+h)-k(x0-h))/(h*2)

#wypisanie przybliżenia wartości pochodnej w pkt x0
print "k'(x0) = ", pochodnaj_Oh, "z dok. O(h)"
print "k'(x0) = ", pochodnaj_Oh2, "z dok. O(h**2)"

#wykreślenie stycznych określonych przez pochodną funkcji
# k z dokładnością O(h) i O(h**2)
pt.figure(4)
wyk = pt.subplot(111)
#Wykreślenie charakterystyki funkcji 1/(sqrt(1 + x^2))
wyk.plot(to, k(to), 'r', label='1/(sqrt(1 + x**2))')
wyk.grid(True)
wyk.hold(True)
#Wykreślenie stycznej do funkcji (x**2)log(x) z dokładnością O(h)
wyk.plot(t, pochodnak_Oh*t + wart_k_x0, 'b', label='Przyb. pochodnej \nfunkcji 1/(sqrt(1 + x**2)) z dok. O(h)')

#Wykreślenie stycznej do funkcji (x**2)log(x) z dokładnością O(h^2)
wyk.plot(t, pochodnak_Oh2 * t + wart_k_x0, 'k', label='Przyb. pochodnej \nfun. 1/(sqrt(1 + x**2)) z dok. O(h^2)')

#dodanie legendy
wyk.legend()
pt.xlabel('Czas [s]', fontsize=12, color='black')
pt.ylabel('Wartosc funkcji')
pt.title('Charakterystyka 1/(sqrt(1 + x**2)) i jej pochodne', fontsize=12, color='black')
pt.show()