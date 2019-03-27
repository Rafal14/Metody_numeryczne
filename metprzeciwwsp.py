#-*-coding: utf-8 -*-


import scipy as sp
import numpy as np
from math import *

print "Rozwiązywanie układu dwóch równań liniowych z dwoma niewiadomymi metodą przeciwnych współczynników"
print "ax + by = c"
print "dx + ey = f\n"

print "Podaj wartości współczyników układu równań\n"

a = float(raw_input("a = "))
b = float(raw_input("b = "))
c = float(raw_input("c = "))
d = float(raw_input("d = "))
e = float(raw_input("e = "))
f = float(raw_input("f = "))

#Sprawdzenie czy pary współczynników a,d oraz b,e są przeciwne
przec_a = a * (-1)
przec_b = b * (-1)

h = c + f



if a==0 and b==0 and d==0 and e==0:
    if c==0 and f==0:
        print "Układ równań jest nieoznaczony"
        print "Istnieje nieskończenie wiele rozwiazań"
    else:
        print "Układ równań jest sprzeczny"
        print "Nie istnieje żadne rozwiązanie"
elif a==d and b==e and c==f:
    print "Układ równań jest nieoznaczony"
    print "Istnieje nieskończenie wiele rozwiazań"
else:
    # Badanie dla układu nieoznaczonego
    k  = abs(a / d)
    k1 = k * d
    k2 = k * e
    k3 = k * f
    if k1==a and k2==b and k3==c:
        print "Układ równań jest nieoznaczony"
        print "Istnieje nieskończenie wiele rozwiazań"
    else:
        if przec_a == d:
            g = b + e
            y = h/g
            x = (f - e*y)/d
            print "Układ równań jest oznaczony"
            print "Rozwiązanie to:\n"
            print "x = ", x
            print "y = ", y
        elif przec_b == e:
            g = a + d
            x = h/g
            y = (f - d*x)/e
            print "Układ równań jest oznaczony"
            print "Rozwiązanie to:\n"
            print "x = ", x
            print "y = ", y
        else:
            ilocz1 = a * d
            ilocz2 = b * e
            bezwg1 = abs(ilocz1)
            bezwg2 = abs(ilocz2)
            #Wykonanie obliczeń względem mniejszego iloczynu
            if bezwg1<bezwg2:
                stala1 = abs(d)
                stala2 = abs(a)
                if ilocz1 > 0:
                    stala1 = stala1 * (-1)
                b1 = b * stala1
                c1 = c * stala1
                e1 = e * stala2
                f1 = f * stala2
                d1 = d * stala2

                m = b1 + e1
                if m!=0 and d1!=0:
                    y=(c1 + f1)/(b1 + e1)
                    x=(f1-e1*y)/d1
                    print "Układ równań jest oznaczony"
                    print "Rozwiązanie to:\n"
                    print "x = ", x
                    print "y = ", y
                else:
                    print "Układ równań jest sprzeczny"
                    print "Nie istnieje żadne rozwiązanie"
            else:
                stala1 = abs(e)
                stala2 = abs(b)
                if ilocz2 > 0:
                    stala1 = stala1 * (-1)
                a1 = a * stala1
                c1 = c * stala1
                d1 = d * stala2
                e1 = e * stala2
                f1 = f * stala2

                mianownik = a1 + d1
                if mianownik != 0 and e1!=0:
                    x  = (c1 + f1)/(a1 + d1)
                    y  = (f1 - d1*x)/e1
                    print "Układ równań jest oznaczony"
                    print "Rozwiązanie to:\n"
                    print "x = ", x
                    print "y = ", y
                else:
                    print "Układ równań jest sprzeczny"
                    print "Nie istnieje żadne rozwiązanie"