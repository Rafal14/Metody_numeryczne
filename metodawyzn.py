#-*-coding: utf-8 -*-


import scipy as sp
import numpy as np
from math import *

print "Rozwiązywanie układu dwóch równań liniowych z dwoma niewiadomymi metodą wyznaczników"
print "ax + by = c"
print "dx + ey = f\n"

print "Podaj wartości współczyników układu równań\n"

a = float(raw_input("a = "))
b = float(raw_input("b = "))
c = float(raw_input("c = "))
d = float(raw_input("d = "))
e = float(raw_input("e = "))
f = float(raw_input("f = "))

#Sprawdź czy można zastosowań metodę wyznczników
if a==0 and b==0 and d==0 and e==0:
    print "Współczynniki a, b, d, e mają wartość zero"
    print "Nie można zastosować metody wyznaczników"
else:
    #Utwórz macierze na podstawie zadanego układu równań
    A  = np.array([[a, b], [d, e]])
    A1 = np.array([[c, b], [f, e]])
    A2 = np.array([[a, c], [d, f]])

    #Oblicz wyznaczniki macierzy A, A1, A2
    detA  = np.linalg.det(A)
    detA1 = np.linalg.det(A1)
    detA2 = np.linalg.det(A2)

    #Określ liczbę rozwiązań na podstawie wartości wyznaczników
    if detA==0:
        if detA1==0 and detA2==0:
            print "Układ równań jest nieoznaczony"
            print "Istnieje nieskończenie wiele rozwiazań"
        else:
            print "Układ równań jest sprzeczny"
            print "Nie istnieje żadne rozwiązanie"
    else:
        x=detA1/detA
        y=detA2/detA
        print "Układ równań jest oznaczony"
        print "Rozwiązaniem układu równań jest para liczb\n",
        print "x = ", x
        print "y = ", y

