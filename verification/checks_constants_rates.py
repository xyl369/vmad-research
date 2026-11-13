import numpy as np
from scipy.special import sph_harm
hbar=1.054571817e-34; c=2.99792458e8; G=6.67430e-11; kB=1.380649e-23
lP=np.sqrt(hbar*G/c**3); tP=lP/c
eta0=c**3/(4*hbar*G); K0=c**4/(8*np.pi*G)
print("lP",lP,"eta0",eta0,"K0",K0, "check K=hbar c eta/2pi",hbar*c*eta0/(2*np.pi))
print("2pi/(hbar c eta0) vs 8piG/c^4", 2*np.pi/(hbar*c*eta0), 8*np.pi*G/c**4)
# DSW rates  Gamma ~ G m^2 d^4 a^5/(hbar c^10) ; dS: a=cH
H0=2.2e-18
for (m,d,lab) in [(1e-14,1e-4,"BMV 1e-14kg,100um"),(1e-3,1e-3,"1g,1mm"),(1.0,1.0,"1kg,1m")]:
    a=c*H0; Gam=G*m**2*d**4*a**5/(hbar*c**10)
    print("DSW dS",lab,Gam)
# Rindler with a=g (formal)
m,d=1e-14,1e-4; g=9.8
print("DSW Rindler a=g (formal)",G*m**2*d**4*g**5/(hbar*c**10))
# dimension check of total-book (5-5): G Q^2/(hbar c^5 R^6) -> units s^-1 * 1/(c^5 R) -> wrong
# Fibonacci fusion dim for n tau anyons on sphere -> F_{n-1}
F=[0,1]
for i in range(30): F.append(F[-1]+F[-2])
phi=(1+5**.5)/2
print("12 tau anyons sphere dim",F[11],"phi^10",phi**10,"log2",np.log2(F[11]),"ln",np.log(F[11]))
# icosahedron area ratio
a=1; AI=5*np.sqrt(3); VI=5*(3+np.sqrt(5))/12
print("A_I/A_sphere",AI/(36*np.pi*VI**2)**(1/3))
# Sakharov/Susskind-Uglum suppression for material scale a
for aa in [1e-10,1e-15,1e-19]:
    print("delta eta/eta ~ (lP/a)^2 for a=",aa,(lP/aa)**2)
# Verlinde-Zurek length fluctuation sqrt(lP L / 4pi)? 
for L in [1.0,4e3]:
    print("VZ dL ~ sqrt(lP L)",L,np.sqrt(lP*L))
# eta-noise AD model: Gamma = (m g d/hbar)^2 sigma^2 tau_c T ; invert bound from SG noise ~1e-11 g/rtHz
m,d,g=1e-14,1e-4,9.8
pref=(m*g*d/hbar)
print("m g d/hbar [rad/s per unit eps]",pref)
S=(1e-11)**2  # eps PSD 1/Hz (one-sided) ~ dg/g
print("Gamma_AD white (pref^2 * S/2)",pref**2*S/2)
# LLR Gdot/G bound vs H0
H0yr=H0*3.156e7
print("H0 per yr",H0yr,"LLR 1e-13/yr in H0 units",1e-13/H0yr)
# Rayleigh check symbolic: equilibrium Lambda Q^2 = 4pi drho R^4 + 8 pi sigma R^3 > 8 pi sigma R^3 (Rayleigh l=2)
# fission n=1: E ~ Q^{4/3}; halves
print("n=1 split ratio 2^(-1/3)",2**(-1/3))
# dS Gibbons-Hawking temperature
TGH=hbar*H0/(2*np.pi*kB); print("T_GH",TGH)
# BB: kg brain dS fluctuation exponent  M c^2/(kB T_GH)
print("BB exponent ~ Mc^2/kT for 1kg",1*c**2/(kB*TGH))
print("S_dS ~ pi c^2/(G hbar H^2)... ", np.pi*c**5/(G*hbar*H0**2))
# TCC lifetime bound (1/H) ln(Mp/H)
Mp_over_H = np.sqrt(hbar*c**5/(8*np.pi*G))/(hbar*H0)
print("TCC T<=", np.log(Mp_over_H)/H0/3.156e7/1e9,"Gyr")
