# -*-coding: utf-8 -*-

# Program do znajdowania pierwiastków równania
# x**3 + x**2 - 3*x - 2 = 0
# w przedziałach [-3, -1] i [-1, 1]
# metodą siecznych

import numpy as np

def f(u):
    '''Oblicza wartość funkcji wielomianu z zadania'''
    wartosc = u ** 3 + u ** 2 - 3 * u - 2
    return wartosc

print "Program do znajdowania pierwiastków równania\n"
print "x**3  +  x**2  - 3*x  - 2 = 0\n"
print "w przedziałach [-3, -1] i [-1, 1]"
print "metodą siecznych."

# wartości ograniczeń przedziału
a = -3.0
b = -1.0
c = 1.0

N = int(raw_input("\nPodaj liczbę kroków "))

#wartość o którą należy zmniejszyć przedział
#gdy wystąpi błąd iteracji spowodowany tym
#że kolejne wartości |f(x)| nie tworzą ciągu
#malejącego
epsilon=0.25

# obliczenie pierwiastka w przedziale [-3, -1]
x2 = a
x1 = b

# utwórz wektor argumentów i wartości funkcji
arg = []
war = []
#sygnalizuje błąd iteracji
flaga =1

for licz in range (0,2):
    arg.append(x2)
    arg.append(x1)

    # oblicznie wartości funkcji wielomianu dla x2 i x1
    wart_x2 = f(x2)
    wart_x1 = f(x1)
    war.append(wart_x2)
    war.append(wart_x1)

    for k in range(0, N):
        #obliczenie wartości bezwzględnych f(x2), f(x1)
        bezwzg_fx2 = abs(wart_x2)
        bezwzg_fx1 = abs(wart_x1)

        #sprawdzenie czy wartości |f(x)| tworzą ciąg rosnący
        #warunek przerywający iterację
        if bezwzg_fx2 < bezwzg_fx1 and k>0:
            flaga=0

        roznica = wart_x1-wart_x2
        #jeżeli ciąg kolejnych wartości |f(x)| jest malejący
        #i różnica wartości f(x2) i f(x1) jest nie 0,
        #i nie uzyskano przybliżenia pierwiastka funkcji
        if flaga!=0 and roznica!=0 and wart_x1!=0:
            #obliczenie wartości x na podstawie poprzednich wartości
            x = x1 - ( (wart_x1 * (x1-x2))/roznica )
            arg.append(x)
            #obliczenie wartości funkcji dla x
            wart_x = f(x)
            war.append(wart_x)
            #zapamiętanie poprzednich wartości
            x2 = x1
            x1 = x
            wart_x2 = f(x2)
            wart_x1 = f(x1)
        #jeżeli wystąpił błąd iteracji spowodowany tym, że
        #ciąg wartości |f(x)| jest niemalejący
        elif flaga == 0:
            arg_w = np.array([arg])
            dl = arg_w.size
            no = dl - 1
            #usuń wartości oprócz dwóch pierwszych elementów
            while no > 1:
                arg.pop(no)
                war.pop(no)
                no = no -1
            #zmniejsz zakres badanego przedziału funkcji
            if arg[0] < 0:
                arg[0] = arg[0] + epsilon
            else:
                arg[0] = arg[0] - epsilon
            if arg[1] < 0:
                arg[1] = arg[1] + epsilon
            else:
                arg[1] = arg[1] - epsilon
            wart_x2=f(arg[0])
            wart_x1=f(arg[1])
            war[0] = wart_x2
            war[1] = wart_x1
            flaga =1
            k=0
    #zmiana granic przedziału dla obliczenia pierwiastka w przedziale [-1,1]
    x2 = b
    x1 = c

    #konwersja do np.array
    wekt_arg = np.array([arg])
    wekt_war = np.array([war])
    #utworzenie macierzy
    macierz = np.concatenate((wekt_arg.T, wekt_war.T), axis=1)
    print " x,    f(x)"
    print macierz
    print "\n"
    l = wekt_arg.size
    for i in range (0,l):
        arg.pop(0)
        war.pop(0)