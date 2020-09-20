import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.colors
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import pi, G, c, hbar, m_n, h
from scipy.interpolate import interp1d

import math

#Constants
M_sun = 1.989 * 10.**(30.) #SI kg
from_mev = 1.602176565 * 10.**(32)  #This is to SI, to cgs its 10^33 
kappa = (8. * pi)**(1./2.)      # equals sqrt(8 pi G) but G = 1


#---------------------------------------------------------

# READ IN THE EOS 

#---------------------------------------------------------

#STT Thesis modified polytrope

data = np.loadtxt("EOS_files/STT_EOS1_SI.dat",skiprows=1,delimiter='\t')
data = data.transpose()

eos_ps = data[0]                            #first column is pressure
eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities

eos_rhos = data[1]                          #second column is denisty
eos_rhog = eos_rhos * G * c**(-2.)          #Use G * c**-2 because this is mass density

# From INGO

"""data = np.loadtxt("EOS_files/EOS_nsat_ind1.dat",skiprows=0,delimiter='\t')
data = data.transpose()
#print np.shape(data)

eos_ps = data[1]      
eos_ps = eos_ps * from_mev                  #second column is pressure
#print max(eos_ps)
eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities

eos_rhos = data[2]    
eos_rhos = eos_rhos * from_mev              #Third column is denisty
eos_rhog = eos_rhos * G * c**(-4.)          #Use G * c**-4 because this is energy density"""

#SLY

data = np.loadtxt("sly_compare/sly_eos_dense.dat",skiprows=1,delimiter='\t')
data = data.transpose()
sly_ps = data[1]                #pressure in cgs
sly_ps = sly_ps * 0.1           #pressure in SI
#print "pressure = ", eos_ps
#eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities


sly_rhos = data[2]              #mass density in cgs
sly_rhos = sly_rhos * 10**3.    #mass 
#print "mass densities = ", eos_rhos
#eos_rhog = eos_rhos * G * c**(-2.)          #Use G * c**-2 because this is mass density




def EOS_geometric_fromfile(p):
    ps = eos_pg
    rhos = eos_rhog

    rho = np.interp(p,ps,rhos)
    return rho

def EOS_scaled_fromfile(p):
    ps = eos_pg/p_c
    rhos = eos_rhog/rho_c

    rho = np.interp(p,ps,rhos)
    return rho





def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    rho = EOS_scaled_fromfile(p)
    
    temp1 = r*(r - 2. * M )
    temp2 = M + 4 * np.pi * r**3.  * p * p_c 

    dMdr = 4. * np.pi * rho * rho_c * r**2.
    dpdr = - 1. * (p * p_c + rho * rho_c) * ( temp2 / temp1 ) * (1./p_c)


    dphidr = temp2 / temp1 



    return [dMdr,dpdr,dphidr]

def STT_TOV(r,y):
    M = y[0]
    p = y[1]
    nu = y[2]
    phi = y[3]
    psi = y[4]
    #Mbar = y[5]    #to get mbar, I need mb and n from the EOS

    
    A = math.exp(alpha_0 * phi)
    A4 = A**(4.)
    kappa2 = kappa**2.
    alpha_tilde =  phi * A 

    rho = EOS_scaled_fromfile(p)
    
    
    temp1 = (r-2.*M)

   
    dMdr = (kappa2/2.) * r**2. * A4 * rho + (1./2.)*r*temp1 * psi**2.

    dp_temp = (kappa2/2.) * ((r**2. * A4 * p)/(temp1)) 
    dp_temp = dp_temp + (1./2.)* r * psi**2.
    dp_temp = dp_temp + M/(r * temp1)
    dp_temp = dp_temp + alpha_tilde*psi
    dpdr = -(rho - p)* dp_temp

    dnudr = kappa2 * ((r**2. * A4 * p)/(temp1)) + r * psi**2. + ((2. * M)/(r*temp1))
    dphidr = psi
    dpsidr = (kappa2/2.) * ((r * A4)/(temp1)) * ( alpha_tilde * (rho - 3.*p) + r * psi * (rho - p))  - ((2*(r-M))/(r*temp1))*psi

    return [dMdr,dpdr,dnudr,dphidr,dpsidr]



def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

alpha_0 = 10.**(-3.)            #Coupling constant
r_0 = 10.**(-5.) #m
r_stop = 40. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

p_low = 1. * from_mev                           #Ingo says to start with ~ 1MeV cm^-3
#print p_low

pressures = np.logspace(30,36,num=2)           #p_max in the file from Ingo is 2.4 x 10^35, so max needs to be less than that 
#pressures  = [10**34.,10**34. 
pressures = pressures * G * c**(-4.)

#print pressures


radii = []
masses = []


nums = np.linspace(10.**5.,10.**8.,num=100)
#step_size = np.linspace(1.,10**(-5.))
step_size = 1.
#print step_size
rho_c_list = []


for x in pressures:

    #Initial conditions
    p_0 = x
    rho_c = EOS_geometric_fromfile(p_0)           #the thesis uses rho_c for both pressure and density scaling
    rho_0 = EOS_geometric_fromfile(p_0)
    p_c = rho_c
    M_0 = 0.
    nu_0 = 0
    phi_0 = 4. * 10**(-3.)
    A4_0 = (math.exp(alpha_0 * phi_0))**4. 
    alpha_til_0 = phi_0 * A4_0
    psi_0 = (4.*pi/3)* r_0 * A4_0 * alpha_til_0 * (rho_c - 3. * p_0)

    rho_c_list.append(rho_c)

    t_eval = np.linspace(r_0,r_stop,num=10.**5.)
    p_0 = x/p_c
    phi_0 = 0.
    #y0 = [M_0,p_0,phi_0] 

    y0 = [M_0, p_0, nu_0, phi_0, psi_0]
    print y0
    #print y0

    #soln = solve_ivp(TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)     #I took out dense_output=True for speed
    soln = solve_ivp(STT_TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)

    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    #print soln

    M = M * c**2. / G
    M = M/ M_sun
    #p = p * p_c



    radii.append(r[-1]/1000.)
    masses.append(M[-1])


    print p[0], r[-1], M[-1]

    """plt.plot(r,M)
    plt.xlabel("radius (m)")
    plt.ylabel("mass")
    plt.show()
    plt.close()
    #plt.savefig("plots/mass_profile_geometric_" + str(i) + ".pdf")
    #plt.close()"""



    #plt.plot(r,p, label='numerical')
    #plt.plot(r,ps,label='calculated')
    #plt.xlim(0,r_stop)
    #plt.ylim(0,x)
    #plt.legend()
    #plt.xlabel("radius (m)")
    #plt.ylabel("pressure")
    #plt.savefig("plots/pressure_profile_SI_" + str(i) + ".pdf")
    #plt.close()
    #i+=1


#f= open("EOS_files/mass_radius_447.dat",'w+')
#f.write("Mass (M_sun)) \t Radius (km) \n")
for x in range(0,len(radii)):
    if radii[x] < 20.0:
        temp = str(masses[x]) + "\t" + str(radii[x]) + "\n"
        #f.write(temp)
        print temp

r = r/1000.


#test cgs
#p = p * p_c #geometric 
#p = p * c**4. / G  # SI
#p = p * 10. #bayre



#print masses 

"""plt.plot(r,p, label='numerical')
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel("pressure")
plt.savefig("plots/GR/pressure_profile_dim.pdf")
plt.close()"""

#print np.shape(radii)
#print np.shape(masses)
"""for x in range(0,len(rho_c_list)):
    rho_c_list[x] = rho_c_list[x] * c**2. / G
    rho_c_list[x] = rho_c_list[x] / 10**(14.)

plt.scatter(rho_c_list,masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())
plt.legend()
#plt.xlim(0,50)
plt.ylim(0,3)
plt.xlabel(r'$\rho_c[10^14 g/cm^3]$')
plt.ylabel(r'Mass (M$_{\odot}$)')
plt.show()"""

plt.plot(sly_ps, sly_rhos, label='sly')
plt.plot(eos_ps, eos_rhos, label='polytrope')
plt.xlabel("pressure")
plt.ylabel('density')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

plt.scatter(masses, radii,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())
plt.legend()
plt.ylabel(r'radius (m)')
plt.xlabel(r'Mass (M$_{\odot}$)')
plt.xlim(0,3)
plt.ylim(0,20)
#plt.savefig("plots/GR/mass_radius_chiraleft_447.pdf")
#plt.close()
#plt.show()

#print r_stop
#print r_0

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


