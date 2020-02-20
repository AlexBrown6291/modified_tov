import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

#Constants
G = 6.67408 * 10.**(-11.) #SI m^3 kg^-1 s^-2
c = 3.00 * 10.**8. #SI m/s
e_0 = 1.603 * 10.**(38.) #cgs   ergs/cm^3 
e_0 = e_0 * 1. * 10.**(-7.) #Convert to J/cm^3
e_0 = e_0 * 100.**(3.) #convert to J/m^3 (SI)
e_0_g = e_0 * G * c**(-4.) # Geometric units
M_sun = 1.989 * 10.**(30.) #SI kg

def EOS(p):
    K = 1.914 
    rho = 1. * p**(-3./5.) / K
    print p
    return rho

def TOV(r,y):
    M = y[0]
    p = y[0]

    rho = EOS(p)

    dMdr = 4 * np.pi * r**2. * rho * e_0_g
    dpdr = - 1. * rho * M  / (r**2.)

    if p<= 0:
        dMdr = 0.
        dpdr = 0.

    return [dMdr,dpdr]


#Initial Conditions

M_0 = 0.
p_0 = 1.  #dimensionless
r_0 = 0.1 #m
r_stop = 10 #km
r_stop = r_stop * 10.**(3.) #SI = m

y0 = [M_0,p_0]
t_span = (r_0,r_stop)

print y0
print t_span

soln = solve_ivp(TOV,t_span,y0,method='RK45')



