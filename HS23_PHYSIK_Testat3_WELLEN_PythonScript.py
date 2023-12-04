
# Saletti - hier wieder ein Lösungsscript für das HS2023 -> es gibt wieder mal einige neue Aufgaben :)
# Ähnliches Prinzip wie beim letzten mal: Die Aufgaben sind in Zellen unterteilt. 
# In den jeweiligen Zellen in die Inputbereiche alles reinhämmern und über ctrl + k die Lösungen 
# im Consolen-Moped rausfetzen lassen 8). Hier und da sind ein paar Randbemerkungen und Notitzen, die evt. beim
# Lösen helfen. 

  'n Gruss und tschau Kakao! - Julian

#%% Aufgabe 1 Rüttler (MEP Aufgabe)
import numpy as np

#-----------------INPUT------------------
f = 44.3#Hz
Fs = 1.2#N
l = 42#cm
#---------------------------------------

#Einheiten umrechnen
l = l*1e-2#m

#A) 
lambda_2 = (2 * l) / 2
c = f * lambda_2
u = (Fs) / (c**2) #1e3 um von kg/m in g/m umzurechnen

print("Aufgabe 1\n A) Längendichte u = ", u*1e3, "g/m")

#B)
lambda_3 = (2 * l) / 3
c = f * lambda_3
Fs = u * c**2

print("B) Spannkraft = ", Fs, "N")


#%% Aufgabe 2 Autbahnraststätte (MEP FS18)
import numpy as np

#-----------------INPUT------------------
L1 = 65#dB
d = 19#m
L2 = 71#dB
#---------------------------------------

I0 = 1e-12#const

#A)
I_tot = 10**(L2 / 10) * I0
I_mitZ = 10**(L1 / 10) * I0

I_Zikade = I_tot - I_mitZ
print("Aufgabe 2\n A) Iz = ", I_Zikade, "W/m^2")

#B)
Pav = I_Zikade * 4 * np.pi * d**2 #Kugelwelle Formel umgestellt
print(" B) Pav = ", Pav, "W")

#C)
#-----------------INPUT------------------
L_alleRastenAus = 89#dB
#---------------------------------------
I_alleRastenAus = 10**(L_alleRastenAus / 10) * I0
n = I_alleRastenAus / I_Zikade
print(" C) n = ", n)

#%% Aufgabe 3 Leben im Tessin
import numpy as np

#-----------------INPUT------------------
d1 = 1.2#km
d2 = 2.9#km
L1 = 66#dB
#---------------------------------------

d1 = d1 * 1e3 #einheiten in m umrechnen
d2 = d2 * 1e3

#A)
L2 = L1 - 20 * np.log10(d2 / d1)
print("Aufgabe 3\n A) L_B = ", L2, "dB")

#B)
#-----------------INPUT------------------
L_mitHund = 68#dB
#---------------------------------------

I0 = 1e-12#W/m^2

I_tot = 10**(L_mitHund / 10) * I0
I_Autobahn = 10**(L1 / 10) * I0

IH = I_tot - I_Autobahn
print(" B) I_H = ", IH* 1e6, "uW / m^2")

#C)
#-----------------INPUT------------------
L_mitZwinger = 74#dB
#---------------------------------------

I_tot = 10**(L_mitZwinger / 10) * I0

IZwinger = I_tot - I_Autobahn

n = IZwinger / IH
print(" C) n = ", n, "Hunds im Zwingsmoped")

#%% Aufgabe 4 Das Krankenauto fährt zu schnell MEP HS 18
import numpy as np

#-----------------INPUT------------------
c = 340#m/s
ue = 50#km/h
f1 = 267#Hz
f2 = 119#Hz
#---------------------------------------

ue = ue * ( 1e3 / 60**2 ) #Einheiten umrechnen in m/s

#A) 
#beide dopplerformeln umstellen nach us

#habe viel substituiert - macht umstellen wesentlich leichter
u = c + ue
v = c - ue
q = ( (f2/f1) * (u / v ) )

us = ( c * (1 - q) ) / ( q + 1 )

f0 = f1 * ( (c - us) / (c + ue) )

print("Aufgabe 3)\n A)\nv = ", us, "m/s\nf = ", f0,"Hz")

#%% Aufgabe 5 Zwei Schallwellen MEP FS 18
import numpy as np

#-----------------INPUT------------------
# Werte aus Funktion ablesen, siehe Harmonische Welle in Formelsammlung

#fn p1:
A1 = 0.014#Pa
k1 = 3.9#1/m
w1 = 1287.0#1/s
phi1 = 1.045#rad

#fn p2
A2 = 0.018#Pa
k2 = 3.85#1/m
w2 = 1270.5#1/s
phi2 = 0
#---------------------------------------

#A) 
vc = w1 / k1 #Schallgeschwindigkeit ist bei beiden gleich, mit welchen faktoren man rechnet also egal
print("Aufgabe 5)\nA) vc = ", vc, "m/s")

#B) 
#-----------------INPUT------------------
Ip1 = 2.5e-7#W/m^2
#---------------------------------------

rho = A1**2 / ( Ip1 * vc * 2 ) #Einheiten gehen in der Form auf
print("B) rho = ", rho, "kg/m^3")

#C)
d_w = np.abs(w1 - w2)
Amax = A1 + A2 #Amplituden addieren sich im Maximum
Amin = np.abs(A2 - A1) # und subtrahieren sich im Minimum

Imax = 0.5 * ( Amax**2 / (rho * vc ) )
Imin = 0.5 * ( Amin**2 / (rho * vc ) )
print("C) Imax = ", Imax, "W/m^2\n   Imin = ", Imin, "W/m^2")
print("darauf achten die Zehnerpotenz mit zu nehmen! :)")

#%% Aufgabe 6 Lautsprechertest im reflexionsfreien Raum MEP HS18
import numpy as np

#-----------------INPUT------------------
c = 340#m/s
d = 16#m
f = 915#Hz
#---------------------------------------

#A)

#-----------------INPUT------------------
LA = 44#dB
#---------------------------------------
I0 = 1e-12#const
# annahme: 2 Kugelquellen ABER nur einer ist aktiv

IA = 10**(LA / 10)*I0
Pav = IA * (d/2)**2 * np.pi * 4
print("Aufgabe 6)\n A) P = ", Pav, "W")

#B)

# die Schalldruckamplitude verdoppelt sich. in Former für Ischall kann man sehen dass wenn sich die
# Amplitude verdoppelt, sich auch die Intensistät verdoppelt. Wir können also direkt die Intensität aus
# Teilaufgabe A verdoppelt und den resultierenden Schallpegel errechnen
L = 10 * np.log10( (2*IA)/I0 ) 
print(" B) L = ", L, "dB")

#C)
#-----------------INPUT------------------
x1 = 7.35#m
x2 = 8.65#m
#---------------------------------------

d_x = x2 - x1

lambda_ = c / f

S1 = 4*np.pi*x1**2
S2 = 4*np.pi*x2**2

I1 = Pav / S1
I2 = Pav / S2

# um ca. 3.5 lambda verschoben! 
# destruktive Interferenz

Ic = I1 - I2

#%% Aufgabe 7 Gespanntes Seil
import numpy as np

#-----------------INPUT------------------
l = 27#m
m = 9#kg
c = 23#m/s
L = 3#m <- Wellenzuglänge
n = 2
P = 1.8#J/s
#---------------------------------------

#A)

u = m / l #längendichte
Fs = c**2 * u
print("Aufgabe 7\n A) Fs = ", Fs, "N")

#B) 

lam = L / n

k = (2*np.pi) / lam
w = c * k
f = (w / (2*np.pi))
w = 2 * np.pi * f

A = np.sqrt( (2 * P) / (u * w**2 * c) )
print("B) A = ", A, "m");

#%% Aufgabe 8 stehende Welle
import numpy as np

#-----------------INPUT------------------
# u(x) = A * sin( -> k <- *x + -> phi <- ) * cos ( -> w <- *t )Werte ablesen
k = 3#1/m
w = 1 / 0.002#rad/s
#---------------------------------------

# A)
f = (w / (2*np.pi))
lambda_ = ( ( 2 * np.pi ) / ( k ) )
v = w / k

print("Aufgabe 8\n A)\n f = ", f, "Hz\n lambda = ", lambda_, "m\n v = ", v, "m/s")

# B)
#-----------------INPUT------------------
L = 349#♦m
#---------------------------------------

n = (2*L) / lambda_ # Faktor 2 da ein Lamba der Wellenlänge einer gesamten periode entspricht -> 2 knöten pro lambda
print("B)\n#Knoten = ", n)

# C)
#-----------------INPUT------------------
x = 22.2#m
#---------------------------------------

next_N = 0;
next_A = 0;

x_modulo = x % (lambda_/2) #Modulo! bricht es runter auf eine halbe Periodenlänge (ein Buckel)

print("Aufgabe C)")
if ( x_modulo > (lambda_ / 4) ): # Modulo ist über die Mitte hinaus, / 4, da lambda 2 knoten enthält! und von denen jeweils will ich die Mitte
    next_N = x + x_modulo
    next_A = x - x_modulo
    print("nach Mitte!")

if ( x_modulo < (lambda_ / 4) ): # Modulo ist noch nicht über die Mitte hinaus! / 4, da lambda 2 knoten enthält!
    next_N = x - x_modulo
    next_A = x + x_modulo
    print("vor Mitte!")

print("Position Wellenbauch: ", next_A, "m\n Position Wellenknoten: ", next_N, "m")

#%% Aufgabe 9 Die Erde bebt MEP 2017
import numpy as np

#-----------------INPUT------------------
t = 0
s_ZA = 10#km
c = 5#km/s
#---------------------------------------

# c = m/s -> Die Welle ist schon 7.5km weit gekommen
# also ist die Zeit errechnebar durch... 

s = 7.5#km abgelesen 

t = s / c
print("Aufgabe 9\n A) t = ", t, "s")

#B) 

#Amplituden verhalten sich antiproportional zu ihrem Abstand zum Epizentrum. Daraus folgt...

#Abgelesen Werte:
A1 = 5 #Amplitude bei x = 4km
d1 = 4#km
d2 = 10#km
lambda_ = 2 #km

A2 = ( d1 / d2 ) * A1
T = lambda_ / c

print("B) \n Ausschlag Ort A: ",A2, "\n Periode = ", T, "s")

#C)

# same wie eben, Intensität nimmt antiproportional zu r^2 zu -> konstante kürzt sich

s_ZB = np.sqrt( s_ZA**2 + 15**2) # pythagoras
Ia_Ib = ( s_ZA**2 / s_ZB**2 )
print("C) Ib / Ib = ", Ia_Ib)

#D)
# beschleunigung -> Wellengleichung 2 mal nach t abgleiten. Trigofunktion kann man ignorieren,
# soll ja um das maximum gehen ( cos = 1 )

# noch ein Paar Einheiten umrechnen
A2 = A2 * 1e-2 #in m
c = c * 1e3 # in m
lambda_ = lambda_ * 1e3 # ebenfalls in m

k = ( 2 * np.pi ) / lambda_
amax = A2 * k**2 * c**2
print("D) a_max = ", amax, "m/s")

#E) 
# die Welle kehrt sich um - man muss also ermitteln, an welchem Wellenabschnitt sie auftrifft
lambda_pos = 10e3 / lambda_
# die welle schwingt volle 5 lambda bis sie an Oberfläche ankommt. 
# ein lambda enthält eine volle Schwingung. der erste Knoten kommt also bei 1/4 lambda

print("E) liegt 0.25 lambda unter A", 0.25 * lambda_ * 1e-3, "km")

#%% Aufgabe 10 Grillieren mit Grillen
import numpy as np

#-----------------INPUT------------------
L_amsel = 60#dB
L_grille = 57#•dB
#---------------------------------------
I0 = 1e-12

#A) wenn man in Intensitäten umrechnet und nachher wieder logarhytmiert kürzt sich recht viel weg 
L_tot = L_amsel * L_grille * 1/10
print("Aufgabe 10\n A) L = ", Ltot, "dB")

#B) 
#-----------------INPUT------------------
L_konzert = 64#dB
#---------------------------------------
I_amsel = I0 * 10**(L_amsel / 10)
I_grille = I0 * 10**(L_grille / 10)
I_konzert = I0 * 10**(L_konzert / 10)

n = (I_konzert - I_amsel) / I_grille

print("B) n = ", n)

#%% Aufgabe 11 Kinderparty mit lauter Musik
import numpy as np

#-----------------INPUT------------------
Pa = 0.3#W
#---------------------------------------
c = 340#m/s - const


#A) 

#-----------------INPUT------------------
da = 5.8#m
db = 4#m
#---------------------------------------

# intensitäten sind gleich 
S1 = 4 * np.pi * da**2
S2 = 4 * np.pi * db**2

Pb = Pa * (S2 / S1)
print("Aufgabe 11\n A) Pb = ", Pb, " W")

#B) 

#-----------------INPUT------------------
F = 250#Hz
#---------------------------------------

# gefragt ist hier die Frequenz, bei der es zu einer destruktiven Imterferenz kommt
# da spassiert bei einer Schiebung von 0.5 lambda

ds = da - db # Streckenunterschied zwischen den Lautsprechern
lambda_ = ds / 0.5 # aus Streckenunterschied benötigtes Lambda rausfinden (diese muss 0.5 lambda entsprechen)
f = c / lambda_
print("B) f = ", f, "Hz")

#C) 

#-----------------INPUT------------------
f = 189#Hz
#---------------------------------------

print("C) lambda: ", c / f) # Hier herrscht absolut konstruktive interferenz! (Lambda ist ganzes vielfaches von ds)
# an Formel für Intensität kann man ablesen, dass Intensität dem Quadrat der Amplitude proportional ist
# Daher is gesamtintensität Summe der Wurzeln der Intensität

Ia = Pa / S1
Ib = Pb / S2

I_tot = Ia + Ib

print(" I = ", I_tot, "W/m^2")
 




#D) 

#-----------------INPUT------------------
#---------------------------------------


























