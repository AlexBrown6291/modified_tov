import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib.colors
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import pi, G, c, hbar, m_n, h
import math

#Constants
"""G = 6.67259 * 10.**(-11.) #SI m^3 kg^-1 s^-2
c = 3.00 * 10.**8. #SI m/s
hbar = 1.0546 * 10.**(-34.) #SI = J * s 
h = hbar*2.*np.pi"""
M_sun = 1.989 * 10.**(30.) #SI kg

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


K_nr = (3.0*pi**2)**(2.0/3.0) * hbar**2. / (5.0*m_h**(8.0/3.0))

K_nr_bar = K_nr * G**(-2./3.) * c**(-2./3.) 

gamma_1 = 5./3.     #non-relativistic degenercy at low densities
gamma_2 = 3.      #relativistic degenercy at high densities

rho_t = 5. * 10.**17.   #Transition Density
rho_t_bar = rho_t * G * c**(-2.)   #to geometric units

p_t = K_nr * rho_t ** gamma_1       #Transition pressure
#p_t_bar = p_t * G * c**(-4.)

p_t_bar = K_nr_bar * rho_t_bar **gamma_1 

K_r_bar = K_nr_bar * rho_t_bar ** (gamma_1 - gamma_2)   #Relativistic K, calculated using continuity



#resulting values of K are the same
#print "c^(-4./3.) = ", c**(-4./3.)
#print "G^(-1/3) = ", G**(-1./3.)

K_bar = K * G**(-1./3.) * c**(-4./3.)

#print "geometric K = ", K_bar



def EOS(p):
    rho = 0.

    if p * p_c < p_t_bar:
        rho = (p_c * p/K_nr_bar)**(1./gamma_1) / rho_c
    else:
        rho = (p_c * p/K_r_bar)**(1./gamma_2) / rho_c


    return rho


def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    rho = EOS(p)
    
    temp1 = r * (r - 2. * M )

    temp2 = M + 4 * np.pi * r**3.  * p * p_c 

    dMdr = 4. * np.pi * rho * rho_c * r**2.
    dpdr = - 1. * (p * p_c + rho * rho_c) * ( temp2 / temp1 ) * (1./p_c)
    dphidr = temp2 / temp1 



    return [dMdr,dpdr,dphidr]



def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

 
M_0 = 0.
r_0 = 0.001 #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

#p = np.logspace(25.,45.,num=10) #SI
pressures = np.logspace(34,45,num=2)
#pressures  = [10**34.,10**34. 
pressures = pressures * G * c**(-4.)

print pressures

i = 0
radii = []
radii2 = []
masses = []
masses2 = []

pressures = np.logspace(20.,40.,num=100)
pressures = pressures * G * c**(-4.)

nums = np.linspace(10.**5.,10.**8.,num=100)
step_size = np.linspace(1.,10**(-5.))
print step_size

for i in (0,1):
    print i
    if i == 1:
        gamma_2 = 2.5 
        K_r_bar = K_nr_bar * rho_t_bar ** (gamma_1 - gamma_2)   #Relativistic K, calculated using continuity
    
    for x in pressures:
        #x = int(x)
        #print x
        p_c = x
        rho_c = (p_c/K_bar)**(1./gamma)

        t_eval = np.linspace(r_0,r_stop,num=10.**5.)
        p_0 = x/p_c
        phi_0 = 0.
        y0 = [M_0,p_0,phi_0]
        #print y0
        #y0 = [x]
        #print y0

        soln = solve_ivp(TOV, t_span, y0,  method='RK45', events=star_boundary, dense_output=True)


        r = soln.t
        M = soln.y[0]
        p = soln.y[1]
        #print soln

        M = M * c**2. / G
        M = M/ M_sun
        #p = p * p_c


        if i == 0:
            radii.append(r[-1]/1000.)
            masses.append(M[-1])
        if i == 1:
            radii2.append(r[-1]/1000.)
            masses2.append(M[-1])

        #print p[0], r[-1], M[-1]

        """plt.plot(r,M)
        plt.xlabel("radius (m)")
        plt.ylabel("pressure")
        plt.savefig("plots/mass_profile_geometric_" + str(i) + ".pdf")
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
        i+=1"""




"""y0 = [0,1.]
p_c = 10.**33. * G * c**(-4.)
rho_c = (p_c/K_bar)**(1./gamma)
soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, t_eval=t_eval, dense_output=True, rtol=10.**(-6.))
r = soln.t
#print soln.t
M = soln.y[0]
p = soln.y[1]
#print p"""

r = r/1000.


#test cgs
#p = p * p_c #geometric 
#p = p * c**4. / G  # SI
#p = p * 10. #bayre


#print masses 

plt.plot(r,p, label='numerical')
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel("pressure")
plt.savefig("plots/GR/pressure_profile_dim.pdf")
plt.close()

plt.scatter(radii,masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm(),label=r'$\gamma=3$')
plt.scatter(radii2,masses2,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm(),label=r'$\gamma=2.5$')
plt.legend()
plt.ylim(0.5,2.5)
plt.xlim(6,18)
plt.xlabel(r'radius (m)')
plt.ylabel(r'Mass (M$_{\odot}$)')
plt.savefig("plots/GR/mass_radius.pdf")
plt.close()

print r_stop
print r_0

"""
plt.scatter(step_size,masses)
plt.xscale('log')
plt.xlim(min(interval)*.5, max(interval)*2)
#plt.ylim(1.42955,1.4296)
plt.legend()
plt.xlabel("Step size (m)")
plt.ylabel("Mass (M_sun)")
plt.savefig("plots/convergence_test.png")
plt.close()"""

"""ps = np.logspace(22.,24,num=200) 
ps = ps * 10.
rhos = EOS(ps)
rhos = rhos / 1000.
rhos = rhos * (3.*10.**10.)**(2.)
ps = ps * 10**(-1.)


plt.plot(rhos,ps)
plt.xlabel("energy density ergs/cm3")
plt.ylabel("pressure dynes")
plt.savefig("plots/pressure-radius_check.pdf")"""


