import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

#Constants
G = 6.67259 * 10.**(-11.) #SI m^3 kg^-1 s^-2
c = 2.99792458 * 10.**8. #SI m/s
hbar = 1.05457266 * 10.**(-34.) #SI = J * s 
h = hbar*2.*np.pi
M_sun = 1.989 * 10.**(30.) #SI kg
#n = 1.
n = 3./2.
gamma = (n+1.)/n
print gamma

m_e = 9.0938215 * 10.**(-31.) # kg
#m_h = 1.6749286* 10.**(-27.) # kg this is neutron mass
m_h =  1.6733 * 10.**(-27.)  #this is the avg hydroigen mass

K = (hbar**2.)/(15 * np.pi **2. * m_e) 
K = K * ((3. * np.pi**2.)/(m_h))**(5./3.)

print K
K  = K * 0.5**(-5./3.)

def EOS(p):

    rho = (p/K)**(1./gamma)

    return rho


def TOV(r,y):
    M = y[0]
    p = y[1]

    rho = EOS(p)
    #print p
    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - G * M * rho / r**2.
    #print dpdr


    return [dMdr,dpdr]

def TOV_odeint(y,r):
    M = y[0]
    p = y[1]

    rho = EOS(p)

    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - rho * G * M * r**(-2.)

    return [dMdr, dpdr]


def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

 
M_0 = 0.
r_0 = 0.01 #m
r_stop = 20000. #km
r_stop = r_stop * 10.**(3.) #SI = m

t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

#p = np.logspace(25.,45.,num=10) #SI
pressures = np.logspace(13,23,num=10)
i = 0
radii = []
masses = []


for x in pressures:
    print x
    y0 = [M_0,x]
    #y0 = [x]
    print y0

    soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, dense_output=True, t_eval=t_eval)


    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    
    radii.append(r[-1])
    masses.append(M[-1])
    print p[-1], r[-1], M[-1]

    plt.plot(r,p)
    plt.xlabel("radius (m)")
    plt.ylabel("pressure")
    plt.savefig("plots/pressure_profile" + str(i) + ".pdf")
    plt.close()
    i = i+1


rho = EOS(pressures)
rs = np.logspace(0,15)
ps = K*rs**gamma

plt.plot(pressures,rho,label='Func')
plt.plot(ps,rs,label='exact')
plt.xlabel("pressure (Pa)")
plt.ylabel("density kg/m^3")
plt.savefig("plots/pressure_density_nr.pdf")
plt.close()


plt.plot(radii,masses)
plt.legend()
plt.xlabel("radius (m)")
plt.ylabel("Mass (M_sun)")
plt.savefig("plots/mass_radius_nr.pdf")
plt.close()


#One single event 
M_0 = 0.
p_0 = 10.**21.
y0 = [M_0,p_0]

r_0 = 0.01 #m
r_stop = 50000. #km
r_stop = r_stop * 10.**(3.) #SI = m

t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,10000)


soln = solve_ivp(TOV,t_span,y0,method='RK45', events=star_boundary, t_eval = t_eval, dense_output=True)

soln_ode = odeint(TOV_odeint,y0,t_eval)
mass = soln_ode[:,0]
press = soln_ode[:,1]
#r = r_eval/100
#r = r/1000

mout = mass[-1]/M_sun
print press[0], mout

r = soln.t
M = soln.y[0]
p = soln.y[1]
r = r/1000.
    
mout = M[-1]/M_sun    
print "output = ", p[0], r[-1], mout

plt.plot(t_eval, press, label='ode_int')
plt.plot(r,p,label='exact')
plt.legend()
plt.xlabel("radius (m)")
plt.ylabel("pressure (Pa)")
plt.savefig("plots/pressure_radius_nr.pdf")
plt.close()