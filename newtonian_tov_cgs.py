import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

g_cgs = 6.67259 * 10.**(-8.)  #cm^3 g^1 s^-1
hbar_cgs = 1.05457266 * 10**(-27.)  #erg s
c_cgs = 3.0 * 10.**10.  # cm /s^1
mh_cgs = 1.6749 * 10.**(-24.) #g 

K_cgs = hbar_cgs*c_cgs
K_cgs = K_cgs /(12*np.pi**2.)
K_cgs = K_cgs * (3 * np.pi **2.)**(4./3.)
K_cgs = K_cgs * (mh_cgs)**(-4./3.)
print K_cgs
#This is correct yay

K_bar = K_cgs * g_cgs**(-1./3.) * c_cgs**(-4./3.)

print K_bar

n = 3./2.
gamma = (n+1)/n

def EOS(p):
    rho = (p/K_cgs)
    rho = rho**(1./gamma)

    return rho


def TOV(r,y):
    M = y[0]
    p = y[1]

    rho = EOS(p)

    dMdr = 4 * np.pi * rho * r**2.
    dpdr = - rho * g_cgs * M * r**(-2.)

    return dMdr, dpdr


M0 = 0
p0 = 10.**25. 
y0 = [M0,p0]

r_0 = 0.1
r_stop = 20. #km
r_stop = r_stop * 10**3. #to m
r_stop = r_stop * 100 #to cm
r = np.linspace(r_0,r_stop,num=4000)


soln = solve_ivp(TOV,(r_0,r_stop),[M0,p0])
rs = soln.t
m = soln.y[0]
p = soln.y[1]

rs = rs / 100
rs = rs / 10**(3.) #back to km

"""soln = odeint(TOV,y0,r)
m = soln[:,0]
p = soln[:,1]
rs = r/100
rs = rs/1000"""

plt.plot(rs,p)
plt.xlabel("radius km")
plt.ylabel("pressure dynes")
plt.savefig("plots/pressure-radius_cgs.pdf")


ps = np.linspace(10.**23,2*10.**24,100)
rhos = EOS(ps)

e_cgs = rhos * c_cgs**2.

plt.plot(e_cgs,ps)
plt.xlabel("energy density ergs/cm3")
plt.ylabel("pressure dynes")
plt.savefig("plots/pressure-density_cgs.pdf")
