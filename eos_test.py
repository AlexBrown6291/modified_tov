import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

#Constants
G = 6.67408 * 10.**(-11.) #SI m^3 kg^-1 s^-2
c = 3.00 * 10.**8. #SI m/s
hbar = 1.0546 * 10.**(-34.) #SI = J * s 
h = hbar*2.*np.pi
M_sun = 1.989 * 10.**(30.) #SI kg
#n = 1.
n = 3.
gamma = (n+1.)/n
print gamma 


m_e = 9.11 * 10.**(-31.) # kg
m_h = 1.67 * 10.**(-27.) # kg

K_si = (hbar * c /(12*np.pi**2.)) * ((3*np.pi)/(m_h))**(4./3.)
K = K_si * G**(-1./3.) * c**(-4./3.)

#K_nr = (1./20.)*(3./np.pi)**(2./3.) * h**2. * (1./m_e) * m_h**(-5./3.)
#K_nr = K_nr * G**(-2./3.) * c**(-2./3.)

#K = K_nr

print "k = ", K

p = np.logspace(35.,45.,num=100)
rho = (p/K_si)**(1./gamma)

p_g = p * G * c**(-4.) # Geometric units
rho_test = (p_g/K)**(1./gamma) #This equation works, that means the issue is with the substitution
print rho_test

rho_g = []

p_0 = p_g[0]
rho_0 = (p_0/K)**(1./gamma)
p_bar = p_g / p_0
rho_bar = ((p_0 * p_bar)/(K))**(1./gamma)/ (rho_0)
rho_g = rho_0 * rho_bar * c**(2.) / G 


"""for x in p_g:
    p_0 = x
    rho_0 = (p_0/K)**(1./gamma)
    print rho_0
    tmp = (((p_0 * x)/(K))**(1./gamma) )/ (rho_0)
    tmp = tmp  * c**(2.) / G 
    rho_g.append(tmp)"""

print rho_g/rho

p_g = p_bar * p_0 * c**4. / G
#rho_g = rho_g * c**(2.) / G 

plt.plot(p,rho,label='SI')
plt.plot(p_g,rho_g,label='geometric')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel('pressure Pascals')
plt.ylabel('density kg/m^3')
plt.savefig('plots/test_eos.pdf')


#Test The Equations in the TOV solver

