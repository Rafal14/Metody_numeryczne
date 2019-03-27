# -*-coding: utf-8 -*-

# Program do znajdowania dodatnich pierwiastków równania
# exp(x) - 3*(x**2) = 0
# metodą Newtona

import numpy as np

def f(u):
    '''Oblicza wartość funkcji z zadania'''
    wartosc = np.exp(u) - 3 * (u**2)
    return wartosc

def pochodnaf(x0, war_h):
    '''Oblicza wartość pochodnej funkcji w zadanym punkcie z dok. O(h**2)'''
    pochf_Oh2 = (f(x0 + war_h) - f(x0 - war_h)) / (war_h * 2)
    return pochf_Oh2


print "Program do znajdowania dodatnich pierwiastków\n"
print "równania exp(x) - 3 * x**2 = 0"
print "metodą Newtona."

# przyjęta wartość punktu końcowego
b = 4.0

h = float(raw_input("\nPodaj wartość h "))

N = int(raw_input("\nPodaj liczbę kroków k "))

#utworzenie wektora argumentów x, wartości f(x) i f'(x)
arg  = []
war  = []
pwar = []
arg.append(b)
war.append(f(b))
pwar.append(pochodnaf(b,h))
#ustal wartość początkową x
x_n = b
#obliczenie
for licz in range(0,2):
    for k in range(0,N):
        #obliczenie wartości x_{n+1} - nowe wartości x
        x_n1 = x_n - ( f(x_n)/pochodnaf(x_n, h) )

        arg.append(x_n1)
        war.append(f(x_n1))
        pwar.append(pochodnaf(x_n1,h))
        #zapamiętaj daną historyczną
        x_n = x_n1

    #utworzenie wektorów kolumnowych
    wek_arg  = np.array([arg])
    wek_war  = np.array([war])
    wek_pwar = np.array([pwar])

    #utworzenie macierzy zwierającej x, f(x) i f'(x)
    macierz = np.concatenate((wek_arg.T, wek_war.T, wek_pwar.T), axis=1)
    print macierz, "\n"

    dlugosc = wek_arg.size
    r = arg[dlugosc-1] / 2
    x_n = r

    #usuń wszystkie wartości z wektorów arg, war, pwar
    for i in range(0, dlugosc):
        arg.pop(0)
        war.pop(0)
        pwar.pop(0)
    #dodaj nową wartości do obliczenia następnego pierwiastka
    arg.append(r)
    war.append(f(r))
    pwar.append(pochodnaf(r, h))

