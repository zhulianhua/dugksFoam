#!/usr/bin/python

from numpy import *

## molecular collision model constant
alpha = 1.0
omega_HS = 0.5
omega_VHS = 0.81

## Boltzmann constant
k = 1.38066e-23

## Argon molecular mass, kg
m = 6.63e-26
## Argon molecular hard sphere diameter, m
d = 4.17e-10
## Argon gas specific gas constant
R = k/m

## Characteristic length, m
L = 1
## Knudsen number
Kn = 0.075


## Characteristic Velocity, m/2
U = 50

## Characteristic temperature, K
T = 273.0

## dimensionless viscosity
muBar = 5.0*(alpha+1.0)*(alpha+2.0)*sqrt(pi)*Kn \
    /(4.0*alpha*(5.0-2*omega_HS)*(7.0-2.0*omega_HS))

## mean free path, m
mfp = L*Kn

## initial density field, Kg/m^3
rho = m/(mfp*sqrt(2.0)*pi*d*d)

## Most probable speed of molecular, m/2
C = sqrt(2*R*T)

## viscosity
mu = rho*C*L*muBar

## Adiabatic constant 
g = 1.4

Ma = U/sqrt(g*R*T)

print "Refer Temperature  (T )  = ", T
print "Knudsen number     (Kn)  = ", Kn
print "Mach number        (Ma)  = ", Ma
print "Reynolds number    (Re)  = ", rho*U*L/mu
print "Viscosity          (mu)  = ", mu
print "Density           (rho)  = ", rho
print "Spefic gas constant (R)  = ", R
print "Most probable speed (C)  = ", C
