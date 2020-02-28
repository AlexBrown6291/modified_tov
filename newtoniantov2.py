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
hbar_g = 137.


"""e_0 = 1.603 * 10.**(38.) #cgs   ergs/cm^3 
e_0 = e_0 * 1. * 10.**(-7.) #Convert to J/cm^3
e_0 = e_0 * 100.**(3.) #convert to J/m^3 (SI)
e_0_g = e_0 * G * c**(-4.) # Geometric units"""

M_sun = 1.989 * 10.**(30.) #SI kg
n = 1.
gamma = (n+1.)/n

m_e = 9.11 * 10.**(-31.)
m_h = 1.67  * 10**(-24.) # g
m_h = m_h / 1000. #kg 
#m_h = m_h * G / c**(2.) #geometric
m_e = m_e * G /c**(2.) #geometric
#print m_e
#K_si = ((hbar*c)/(12 * np.pi)) * ((3*np.pi**2.)/(m_h * c**2.))**(4./3.)
K = (137./(12*np.pi**2.))* ((3*np.pi)/(m_h))**(4./3.)
#print "electron mass ^ 4/3 ", (m_e)**(4./3.)
#print "nucleon mass = ", m_h
#print "K from paper = ", (137./(12*np.pi**2.))* ((3*np.pi)/(m_h))**(4./3.)

K = 20000.
#K = K_si
def EOS(p):
    #K = 1.914 

    rho = ((p_c * p)/(K))**(1./gamma) / rho_c 
    #print p,rho 
    return rho

def TOV(r,y):
    M = y[0]
    p = y[1]

    rho = EOS(p)
    #print rho

    dMdr = 4 * np.pi * r**2. * rho * rho_c 
    dpdr = - ( 1. * rho * M * rho_c) / (r**2. * p_c)

    if p<= 0:
        dMdr = 0.
        dpdr = 0.

    return [dMdr,dpdr]

def star_boundary(r,y):
    return y[1]

star_boundary.terminal = True
#Initial Conditions

M_0 = 0.
#p_0 = 5.  #dimensionless
r_0 = 0.1 #m
r_stop = 20 #km
r_stop = r_stop * 10.**(3.) #SI = m

radii = []
masses = []
t_span = (r_0,r_stop)
#print t_span

p = np.logspace(40.,45.,num=10) #SI
print p
p = p * G * c**(-4.) # Geometric units
print p


i = 0
for x in p:
    p_c = x
    p_0 = x/p_c
    rho_c = (p_c/K)**(1./gamma)

    y0 = [M_0,p_0]
    print "rho_scale = ", rho_c
    print "pscale = ", p_c

    soln = solve_ivp(TOV,t_span,y0,method='RK45',events=star_boundary,dense_output=True)

    r = soln.t
    M = soln.y[0]
    p = soln.y[1]

    M  = M * c**2./G #from geom to SI
    M =  M / M_sun #to solar masses
    #r = r * 10**(-3.) #to km

    z = soln.t_events
    #print z

    radii.append(r[-1])
    masses.append(M[-1])

    """plt.plot(r,M)
    plt.xlabel("radius (km )")
    plt.ylabel("Mass (M_sun)")
    plt.savefig("plots/test_"+ str(x) + ".pdf")
    plt.close()"""


    #Calculate pressure radius curve

    ps = []
    r_n = (((n+1.) * K )/(1.* 4.*np.pi)) * rho_c**(1./n-1.)  #G= 1
    r_n = r_n**(1./2.)
    #print r_n
    #print r
    xi = r / r_n
    print "xi = ", max(xi)
    for j in xi:
        ps.append((math.sin(j)/j)**(n+1.))


    plt.plot(r,p, label='numerical')
    plt.plot(r,ps,label='calculated')
    plt.legend()
    plt.xlabel("radius")
    plt.ylabel("pressure")
    plt.savefig("plots/pressure_profile_pc_" + str(i) + ".pdf")
    plt.close()
    i+=1



#plt.plot(rs,ma, label="func")
plt.plot(radii,masses,label="derived")
plt.legend()
plt.xlabel("radius (m )")
plt.ylabel("Mass (M_sun)")
plt.savefig("plots/mass_radius.pdf")

pressures = []
rhos = []
for x in range(1,20):
    pressures.append(x)
    rhos.append(EOS(x))


"""p = []
rh = []

rh = np.linspace(0,1,20)
for x in rh:
    p.append(1.914 * x**(5./3.))"""

#K = 1.914 
#R = ((np.pi*K)/(2.))**(1./2.)
#R = R*10.**(-3.)
print "exact solution gives r = ", R
print "code solution gives r = ", radii[0]

plt.close()
"""plt.plot(rhos,pressures, label='eos')
#plt.plot(rh,p,label="exact")
plt.legend()
plt.xlabel("density")
plt.ylabel("pressure")
plt.savefig("plots/pressure_density.pdf")"""

