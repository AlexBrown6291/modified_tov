import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math
np.seterr(invalid='ignore')

G = 6.67 * 10.**(-11.) # m^3/kg s^2
c = 3.0 * 10.**8 # m/s

e_0 = 1.603 * 10. ** 38. #cgs
e_0 = e_0 * 10.**-7. #convert to J/cm^3
e_0 = e_0 * 100.**3. #convert to J/m^3
e_0 = G * c**(-4.) * e_0 #Convert to G = c = 1 units

M_sun = 1.989 * 10.**(30.) # kg, SI

def EOS(p):
    K = 1.914 
    n = 1.5
    rho = 1. * (p/K)**(n/(n+1.)) * e_0**(-1./(n+1))
    return rho

def TOV(r,z):
    M = z[0]
    p = z[1]
    
    rho = EOS(p)

    #print "pressure = ", p
    dMdr =  4 * np.pi * (r**2.) * rho * e_0 
    dpdr = - (1. * rho * M )/(r**2.)
    #print "r = ", r
    #print "pressure = ", p
    #print "mass = ", M

    if p<=0:
        dMdr=0
        dpdr=0
        
    
    return [dMdr,dpdr]


#Add stopping condition

def star_boundary(r,y):
    return y[1]

star_boundary.terminal = True




#Initials
r_step =  0.01 #1.0*10**(-5.) # SI m
M_0 = 0. # G = c = 1
p_0 = 20. # unitless = p(Geom) * e_0
r_0 = r_step

#ic = M_0, p_0

r_stop = 20. * 10.**3. 
n = int((r_stop - r_0)/r_step)
print n
r_span = r_0,r_stop
print r_span

pressures = np.linspace(1.,100.,100)
print pressures
#pressures = [1,2]
Masses = []
Radii = []
for x in pressures:
    ic = M_0, x
    soln =  solve_ivp(TOV, r_span, ic, method='RK45',events=star_boundary, dense_output=True)

    t = soln.t
    soln = soln.y
    M = soln[0]
    p = soln[1]
    #print p
    rho = []

    M = M * c**2. / G #to SI=
    M = M/ M_sun
    #print max(M)
    Masses.append(M[-1])

    t = t * 10**(-3.) # radius in km

    Radii.append(t[-1])

    """plt.plot(t,p,label='pressure')
    plt.xlabel('radius (km)')
    plt.ylabel('pressure 1/epsilon_0')
    #plt.ylim(0,1.2)
    plt.savefig("plots/pressure_" + str(x) + ".pdf")
    plt.close()"""



#plt.plot(t,rho,label='density')
"""plt.plot(t,p,label='pressure')
plt.xlabel('radius (km)')
plt.ylabel('pressure 1/epsilon_0')
#plt.ylim(0,1.2)
plt.savefig('plots/pressure.pdf')
plt.close()"""
plt.plot(Radii,Masses)
plt.xlabel('radius (km)')
plt.ylabel('Mass (Solar masses)')
plt.legend()
plt.savefig('plots/with_cutoff_massradius.pdf')