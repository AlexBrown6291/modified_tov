import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math
np.seterr(invalid='ignore')

def EOS(P):
    K = 1.914
    c = 1.
    rho =  (c**2.) * P **(-3./5.)/K
    
    return rho

def TOV(r,z):
    M = z[0]
    p = z[1]
    G = 6.77 * 10.**(-11.) # m^2/kg s
    c = 3.0 * 10.**8 # m/s
    e_0 = 1.603 * 10. ** 38. #SI
    e_0 = e_0 * 10.**-7. #convert to J/cm^3
    e_0 = e_0 * 100.**3. #conver to J/m^3
    #print "e_0 in SI = ", e_0
    e_0 = G * c**(-4.) * e_0 #Convert to G = c = 1 units
    
    rho = EOS(p)

    #print "pressure = ", p
    dMdr = 4 * np.pi * (r**2.) * rho * e_0 
    dpdr = - (1. * rho * M )/(r**2.)

    return [dMdr,dpdr]

#constants
G = 6.77 * 10.**(-11.) # m^2/kg^2
c = 3.0 * 10.**8 # m/s
e_0 = 1.603 * 10. ** 38. #SI
M_sun = 1.989 * 10.**(30.) # kg, SI


#Initials
r_step =  0.01 #1.0*10**(-5.)
M_0 = 0.
p_0 = 10.
r_0 = r_step

ic = M_0, p_0

r_stop = 12. * 10.**3. 
n = int((r_stop - r_0)/r_step)
print n
r = np.linspace(r_0,r_stop,n)
print max(r)
r_span = r_0,r_stop

soln =  solve_ivp(TOV, r_span, ic, method='RK45',t_eval=r, dense_output=True)
#soln = RK45(TOV,r_0,ic,12,max_step=r_step,dense_output=True)

#print soln
t = soln.t
print max(t)
soln = soln.y
M = soln[0]
p = soln[1]
print p
rho = []

for x in p:
    #print x
    rho.append(EOS(x))
    #print EOS(x)

print "MASS ", M
M = M * c**2 / G #to SI=
M = M/ M_sun

t = t * 10**(-3.) # radius in km
#p = p * e_0 #SI
print p


#plt.plot(t,rho,label='density')
plt.plot(t,p,label='pressure')
plt.xlabel('radius (km)')
plt.ylabel('pressure 1/epsilon_0')
#plt.ylim(0,1.2)
plt.savefig('plots/pressure.pdf')
plt.close()
plt.plot(t,M)
plt.xlabel('radius (km)')
plt.ylabel('Mass')
plt.legend()
plt.savefig('plots/newtonian_mr.pdf')