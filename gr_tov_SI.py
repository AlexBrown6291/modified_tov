import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import pi, G, c, hbar, m_n, h
import math


print pi, G, c, hbar, m_n, h 

M_sun = 1.98892 * 10.**(30.) #SI kg
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
K_nr = (3.0*pi**2)**(2.0/3.0) * hbar**2. / (5.0*m_h**(8.0/3.0))

gamma_1 = 5./3.     #non-relativistic degenercy at low densities
gamma_2 = 3.      #relativistic degenercy at high densities

rho_t = 5. * 10.**17.   #Transition Density
p_t = K_nr * rho_t ** gamma_1       #Transition pressure
K_r = p_t / rho_t ** gamma_2     #Relativistic K, calculated using continuity



#resulting values of K are the same
#print "c^(-4./3.) = ", c**(-4./3.)
#print "G^(-1/3) = ", G**(-1./3.)

K_bar = K * G**(-1./3.) * c**(-4./3.)
print "dimensionless K = ", K_bar



def EOS(p):
    rho = 0.
    if p < p_t:
        rho = (p/K_nr)**(1./gamma_1)
    else:
        rho = (p/K_r)**(1./gamma_2)

    return rho


def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    temp = r - 2. * G * M * c**(-2.)
    temp = temp * r

    temp2 = 4. * np.pi * G * r**3. * p * c**(-2.) + G * M 

    rho = EOS(p)
    #print p
    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - (rho + p * c**(-2.)) * temp2 / temp
    dphidr = temp2 / temp 
    #print dpdr


    return [dMdr,dpdr,dphidr]



def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

phi_0 = 0. 
M_0 = 0.
r_0 = 0.01 #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

#p = np.logspace(25.,45.,num=10) #SI
pressures = np.logspace(30,45,num=40)
i = 0
radii = []
masses = []
for x in pressures:
    y0 = [M_0,x,phi_0]
    #y0 = [x]
    print y0
    print "rho = ", EOS(x)

    soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, dense_output=True,max_step=1.)


    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    phi = soln.y[2]
    #print soln

    radii.append(r[-1])
    masses.append(M[-1])

    print p[0], r[-1], M[-1]/M_sun
    print len(M)

    """plt.plot(r,M)
    plt.xlabel("radius (m)")
    plt.ylabel("pressure")
    plt.savefig("plots/GR/mass_profile_SI_" + str(i) + ".pdf")
    plt.close()



    plt.plot(r,p)
    plt.xlim(0,r_stop)
    plt.xlabel("radius (m)")
    plt.ylabel("pressure")
    plt.savefig("plots/GR/pressure_profile_SI_" + str(i) + ".pdf")
    plt.close()"""
    i+=1


"""y0 = [0,10.**33.]
soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, t_eval=t_eval, dense_output=True)
r = soln.t
#print soln.t
M = soln.y[0]
p = soln.y[1]
#print p"""

r = r/1000.


#test cgs

p = p * 10. #bayre


plt.plot(r,p, label='numerical')
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel("pressure")
plt.savefig("plots/GR/pressure_profile_SI.pdf")
plt.close()

for x in range(0,len(masses)):
    masses[x] = masses[x]/M_sun

#radii = radii * 1000.
plt.scatter(radii,masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())
plt.legend()
plt.ylim(0.5,2.5)
plt.xlabel("radius (m)")
plt.ylabel("Mass (M_sun)")
plt.savefig("plots/GR/mass_radius_SI.pdf")
plt.close()



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
