import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

#Constants
G = 6.67259 * 10.**(-11.) #SI m^3 kg^-1 s^-2
c = 3.00 * 10.**8. #SI m/s
hbar = 1.0546 * 10.**(-34.) #SI = J * s 
h = hbar*2.*np.pi
M_sun = 1.989 * 10.**(30.) #SI kg
#n = 1.
n = 3.
gamma = (n+1.)/n
#print gamma

m_e = 9.1094 * 10.**(-31.) # kg
m_h = 1.6749 * 10.**(-27.) # kg

K = ((hbar * c) /(12*np.pi**2.)) * ((3*np.pi**2.)/(m_h))**(4./3.)  #from sanjay's paper
#print K
K = K * 0.5**(4./3.)

#from princeton lecture
K = (h * c)/ 8. 
K = K * (3/np.pi)**(1./3.)
K = K * (1./(m_h))**(4./3.)
#print K
K = K * 0.5**(4./3.)
#print K



#resulting values of K are the same
#print "c^(-4./3.) = ", c**(-4./3.)
#print "G^(-1/3) = ", G**(-1./3.)

K_bar = K * G**(-1./3.) * c**(-4./3.)
print "dimensionless K = ", K_bar



def EOS(p):

    rho = (p/K)**(1./gamma)

    return rho


def TOV(r,y):
    M = y[0]
    p = y[1]

    rho = EOS(p)
    #print p
    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - G * M * rho /r**2.
    #print dpdr


    return [dMdr,dpdr]

"""def TOV(r,y):
    p = y

    #rho = p
    #print p
    #dpdr = - G * np.pi * rho * r**(-2.)
    dpdr = -.5 * p

    #print dpdr

    return dpdr"""

def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

 
M_0 = 0.
r_0 = 0.1 #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

#p = np.logspace(25.,45.,num=10) #SI
pressures = np.logspace(30,40,num=10)
i = 0
radii = []
masses = []
for x in pressures:
    y0 = [M_0,x]
    #y0 = [x]
    print y0

    soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, dense_output=True)


    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    #print soln

    radii.append(r[-1])
    masses.append(M[-1])

    print p[0], r[-1], M[-1]/M_sun

    plt.plot(r,M)
    plt.xlabel("radius (m)")
    plt.ylabel("pressure")
    plt.savefig("plots/mass_profile_SI_" + str(i) + ".pdf")
    plt.close()



    #plt.plot(r,p, label='numerical')
    #plt.plot(r,ps,label='calculated')
    #plt.xlim(0,r_stop)
    #plt.ylim(0,x)
    #plt.legend()
    #plt.xlabel("radius (m)")
    #plt.ylabel("pressure")
    #plt.savefig("plots/pressure_profile_SI_" + str(i) + ".pdf")
    #plt.close()
    i+=1


y0 = [0,10.**33.]
soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, t_eval=t_eval, dense_output=True)
r = soln.t
#print soln.t
M = soln.y[0]
p = soln.y[1]
#print p

r = r/1000.


#test cgs

p = p * 10. #bayre


plt.plot(r,p, label='numerical')
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel("pressure")
plt.savefig("plots/pressure_profile_SI.pdf")
plt.close()


"""radii = radii * 1000.
plt.plot(radii,masses)
plt.legend()
plt.xlabel("radius (m)")
plt.ylabel("Mass (M_sun)")
plt.savefig("plots/mass_radius_SI.pdf")
plt.close()

ps = np.logspace(22.,24,num=200) 
ps = ps * 10.
rhos = EOS(ps)
rhos = rhos / 1000.
rhos = rhos * (3.*10.**10.)**(2.)
ps = ps * 10**(-1.)


plt.plot(rhos,ps)
plt.xlabel("energy density ergs/cm3")
plt.ylabel("pressure dynes")
plt.savefig("plots/pressure-radius_check.pdf")"""




"""p_cgs = np.logspace(23,25,num=100)
e_cgs = p_cgs/K_cgs
e_cgs = e_cgs**(1/gamma)"""

"""e_cgs = np.linspace(10.**27.,20.*10.**27.,200)
rho_cgs = e_cgs / c_cgs**2.
p_cgs = K_cgs * rho_cgs**gamma

plt.plot(e_cgs,p_cgs)
plt.xlabel("energy density ergs/cm3")
plt.ylabel("pressure dynes")
plt.savefig("plots/pressure-radius_check.pdf")"""
