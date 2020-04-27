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

n = 3.
gamma = (n+1)/n
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
r_stop = 25. #km
r_stop = r_stop * 10.**3. #to m
r_stop = r_stop * 100. #to cm
r_eval = np.linspace(r_0,r_stop,num=4000)

radii = []
masses = []

pressures = np.logspace(34,45,num=11) #Use for relativsitic case
#pressures = np.logspace(12,22,num=10) #Use for nonrelativsitic case
r_start = np.logspace(0.00001,0.1,20)
i = 1
#pressures = [10.**34.]
for x in r_start:
    M0 = 0
    p0 = 10.**37.
    y0 = [M0,p0]
    print p0
    soln = solve_ivp(TOV, (x,r_stop), [M0,p0], events=star_boundary, dense_output=True)
    rs = soln.t
    m = soln.y[0]
    p = soln.y[1]

    rs = rs / 100.
    rs = rs / 10**(3.) #back to km

    """soln = odeint(TOV_odeint,y0,r_eval)
    mass = soln[:,0]
    press = soln[:,1]
    r = r_eval/100
    r = r/1000"""

    print rs[-1]
    radii.append(rs[-1])
    masses.append(m[-1])
    """#plt.plot(r,press,label="ode_int")
    plt.plot(rs,p) #,label="solve_ivp")
    plt.legend()
    plt.xlabel("radius km")
    plt.ylabel("pressure bayre")
    plt.savefig("plots/pressure-radius_cgs" + str(i) +".pdf")
    plt.close()"""
    i += 1

for x in range(0,len(masses)):
    masses[x] = masses[x]/m_sun_cgs
    print masses[x]

plt.scatter(radii,masses)
plt.xlabel("radius km")
plt.ylabel("mass g")
plt.savefig("plots/mass-radius_rstart.pdf")
plt.close()





"""plt.plot(rs,m)
plt.legend()
plt.xlabel("radius km")
plt.ylabel("Mass g")
plt.savefig("plots/mass-radius_cgs.pdf")
plt.close()"""



ps = np.logspace(10,41,num=100)
#ps = np.linspace(10.**23,20.*10.**23.,num=100)
#print ps
rhos = EOS(ps)

#e_cgs = rhos * c_cgs**2.

"""plt.loglog(rhos,ps)
#plt.xscale('log')
plt.xlabel("density g/cm3")
#plt.yscale('log')
plt.ylabel("pressure dynes/cm^2")
plt.savefig("plots/pressure-density_cgs_notes.pdf")
plt.close()"""

"""e = np.linspace(10.**27.,20*10.**27,num=100)
e = e/(c_cgs**2.)
p = K_cgs * e**(gamma)
e = e * c_cgs**2.

plt.plot(e,p)
plt.xlabel("energy density ergs/cm3")
plt.ylabel("pressure dynes/cm^2")
plt.savefig("plots/pressure-energydensity_cgs_2.pdf")"""

#One single event 
M_0 = 0.
p_0 = 10.**34.
y0 = [M_0,p_0]

r_0 = 0.01 #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
r_stop = r_stop * 100 #cm 

t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,10000)


soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, dense_output=True, t_eval=t_eval)

r = soln.t
M = soln.y[0]
p = soln.y[1]
r = r/100.
r = r/1000.
    
plt.plot(r,p)
#plt.legend()
plt.xlabel("radius (km)")
plt.ylabel("pressure")
plt.savefig("plots/pressure_profile_cgs.pdf")
plt.close()

mout = M[-1]/m_sun_cgs    
print "output = ", p[0], r[-1], mout