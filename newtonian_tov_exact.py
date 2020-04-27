import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

#Constants 

g_cgs = 6.67259 * 10.**(-8.)  #cm^3 g^1 s^-1
hbar_cgs = 1.05457266* 10**(-27.)  #erg s
h_cgs = hbar_cgs*2.*np.pi
c_cgs = 2.99792458* 10.**10.  # cm /s^1
mh_cgs = 1.6749 * 10.**(-24.) #g 
me_cgs = 9.1093897 * 10.**(-28.) # g
m_sun_cgs = 1.989 * 10.**(33.) #g

K_cgs = hbar_cgs*c_cgs
K_cgs = K_cgs /(12*np.pi**2.)
K_cgs = K_cgs * (3 * np.pi **2.)**(4./3.)
K_cgs = K_cgs * (mh_cgs)**(-4./3.)
print K_cgs
# z/a = 0.5
K_cgs = K_cgs * 0.5 **(4./3.)
print K_cgs
#This is correct yay

#from princeton lecture
K = (h_cgs * c_cgs)/ 8. 
K = K * (3/np.pi)**(1./3.)
K = K * (1./(mh_cgs))**(4./3.)
print K
K = K * 0.5**(4./3.)
#print K

#for gamma = 5/3 
K_nr = hbar_cgs**2.
K_nr = K_nr /(15. *np.pi**2. * me_cgs) 
K_nr = K_nr * (3. * np.pi**2. )**(5./3.)
K_nr = K_nr * (mh_cgs)**(-5./3.)
print "K_non rel = ", K_nr
K_nr = 0.991 * 10.**13.
print "K_non rel = ", K_nr
K_nr = K_nr * 0.5**(5./3.)
#print "K_non rel = ", K_nr
#print "difference = ", 0.991 * 10.**13. - K_nr


"""K_bar = K_cgs * g_cgs**(-1./3.) * c_cgs**(-4./3.)
print "c^(-4/3) = ", c_cgs**(-4./3.)
print "G^(-1/3) = ", g_cgs**(-1./3.)
print "dimensionless k = ", K_bar"""

n = 1.
#n = 5.
gamma = (n+1)/n
K = K_nr
#gamma = 5./3.

def EOS(p):
    rho = (p/K)
    rho = rho**(1./gamma)

    return rho

def EOS_nr(p):
    gamma = 5./3.
    rho = (p/K_nr)
    rho = rho**(1./gamma)

    return rho

def TOV(r,y):
    M = y[0]
    p = y[1]

    #rho = EOS(p)
    rho = EOS(p)

    dMdr = 4 * np.pi * rho * r**2.
    dpdr = - rho * g_cgs * M * r**(-2.)

    return [dMdr, dpdr]

def TOV_odeint(y,r):
    M = y[0]
    p = y[1]

    rho = EOS(p)

    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - rho * g_cgs * M * r**(-2.)

    return [dMdr, dpdr]

def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True


M0 = 0.
p0 = 10.**35. 
y0 = [M0,p0]

r_0 = 0.1
r_stop = 90000. #km
r_stop = r_stop * 10.**3. #to m
r_stop = r_stop * 100. #to cm
r_eval = np.linspace(r_0,r_stop,num=4000)

t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,10000)

soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, dense_output=True, t_eval=t_eval)

r = soln.t
M = soln.y[0]
p = soln.y[1]
    
mout = M[-1]/m_sun_cgs    
print "output = ", p[0], r[-1], 

#Exact Solution for the Radius for n=1

rho0 = EOS(p0)
R = (((1.+1.)*p0 )/(4*np.pi * g_cgs * rho0**2))**0.5 * np.pi  #n=1
#R = (((1.+1.)*p0 )/(4*np.pi * g_cgs * rho0**2))**0.5 * 6.**0.5 #n=0
print "Exact Radius = ", R
#Exact solution of the pressure
print "R_calc - R_exact = ", r[-1] - R
error = abs((r[-1] - R)/R)
print "percent error is ", error

R_s = (((1.+1.)*p0 )/(4*np.pi * g_cgs * rho0**2))**0.5
xi = r / R_s
print max(xi)
theta = np.sin(xi) / xi
P = p0 * theta**(n+1)

#Exact Solution for the Radius for n=5
"""rho0 = EOS(p0)
R_s = (((1.+1.)*p0 )/(4*np.pi * g_cgs * rho0**2))**0.5 
xi = r/R_s
theta = (1.+((xi**2.)/3.))**(-0.5)
P = p0 * theta**(n+1)"""

r = r/100.
r = r/1000.


#plt.plot(r,p,label='numerical')
plt.plot(r,P,label='exact')

plt.xlabel("Radius (km)")
plt.ylabel("pressure dynes/cm^2")
plt.legend()
plt.savefig("plots/pressure-radius_exact.pdf")
plt.close()