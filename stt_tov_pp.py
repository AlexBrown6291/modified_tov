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

"""data = np.loadtxt("EOS_files/STT_EOS1_SI.dat",skiprows=1,delimiter='\t')
data = data.transpose()

eos_ps = data[0]                            #first column is pressure
eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities

eos_rhos = data[1]                          #second column is denisty
eos_rhog = eos_rhos * G * c**(-2.)          #Use G * c**-2 because this is mass density"""

#plt.plot(eos_pg,eos_rhog)



# From INGO

data = np.loadtxt("EOS_files/EOS_nsat_ind1.dat",skiprows=0,delimiter='\t')
data = data.transpose()
#print np.shape(data)

eos_ps = data[1]      
eos_ps = eos_ps * from_mev                  #second column is pressure
#print max(eos_ps)
eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities

eos_rhos = data[2]    
eos_rhos = eos_rhos * from_mev              #Third column is denisty
eos_rhog = eos_rhos * G * c**(-4.)          #Use G * c**-4 because this is energy density

plt.plot(eos_pg,eos_rhog)

#SLY

"""data = np.loadtxt("sly_compare/sly_eos_dense.dat",skiprows=1,delimiter='\t')
data = data.transpose()
sly_ps = data[1]                #pressure in cgs
sly_ps = sly_ps * 0.1           #pressure in SI
#print "pressure = ", eos_ps
eos_pg = sly_ps * G * c**(-4.)              #conversion to dimensionless quantities


sly_rhos = data[2]              #mass density in cgs
sly_rhos = sly_rhos * 10**3.    #mass 
#print "mass densities = ", eos_rhos
eos_rhog = sly_rhos * G * c**(-2.)          #Use G * c**-2 because this is mass density


plt.plot(eos_pg,eos_rhog)"""



#---------------------------------------------------------

# The Bizarre EOS from the paper

#---------------------------------------------------------
#STT Values
m_b = 1.66 * 10.**(-27)  #kg
n_0 = 0.1   #fm^-1

gamma_1 = 2.00
gamma_2 = 2.34
gamma_3 = 2.46
K1 = 0.1 #/ c**2.
K2 = 0.0195 #/ c**2.
K3 = 0.00936 #/ c**2. 

def EOS_pp_1(p):
    #A_num = m_b * n_0
    A_tmp = (K1 * m_b * n_0)**(1./gamma_1)
    A_tmp = A_tmp * (rho_c)**(1.-(1./gamma_1))
    A_eos = m_b * n_0 / A_tmp

    rho = p/(gamma_1 - 1.) * A_eos * p**(1/gamma_1)

    return rho

def get_central_pressure(rho):


#---------------------------------------------------------

# Define Functions

#---------------------------------------------------------



# Equation of State in G=c=1 units
def EOS_geometric_fromfile(p):
    ps = eos_pg
    rhos = eos_rhog

    rho = np.interp(p,ps,rhos)
    return rho

# Equation of State in G=c=1 units with scaled variables
def EOS_scaled_fromfile(p):
    ps = eos_pg/p_0
    rhos = eos_rhog/rho_0

    rho = np.interp(p,ps,rhos)
    return rho

# TOV Equations in GR
def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    rho = EOS_scaled_fromfile(p)
    
    temp1 = r*(r - 2. * M )
    temp2 = M + 4 * np.pi * r**3.  * p * p_0 

    dMdr = 4. * np.pi * rho * rho_0 * r**2.
    dpdr = - 1. * (p * p_0 + rho * rho_0) * ( temp2 / temp1 ) * (1./p_0)


    dphidr = temp2 / temp1 



    return [dMdr,dpdr,dphidr]

# TOV Equations in STT
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
    alpha_tilde =  alpha_0

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
r_start = 10.**(-5.) #m
r_stop = 40. #km
r_stop = r_stop * 10.**(3.) #SI = m

p_low = 1. * from_mev                           #Ingo says to start with ~ 1MeV cm^-3
#print p_low

#pressures = np.logspace(30,36,num=2)           #p_max in the file from Ingo is 2.4 x 10^35, so max needs to be less than that 
pressures = np.logspace(33.5,36,num=50)
#pressures  = [10**34.,10**34. 
pressures = pressures * G * c**(-4.)

#print pressures


radii = []
masses = []


nums = np.linspace(10.**5.,10.**8.,num=100)
#step_size = np.linspace(1.,10**(-5.))
step_size = 0.1
#print step_size
rho_c_list = []

#----- Plot central pressure and densities 

for x in pressures:
    rho_c_list.append(EOS_geometric_fromfile(x))

"""plt.plot(pressures,rho_c_list)
plt.xlabel("pressure")
plt.ylabel('density')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()"""

# variable_c  indicates the central value
# variable_0  indicates a scaling factor


for x in pressures:

    #Initial conditions
    p_c = x                                       # The central density (In geometric units)
    print "p_c = " + str(x) + "in g=c=1 units"
    rho_c = EOS_geometric_fromfile(p_c)           # The central presssure
    print "rho_c = " + str(rho_c) +  "in g=c=1 units"
    rho_0 = rho_c                                 # The scale factor of density 
    p_0 = rho_0                                   # the thesis uses rho_c for both pressure and density scaling
    
    p_c = x/p_0                                   # I need to pass the code scaled values p = p_0 p hat      

    print "p_c = " + str(p_c) + " scaled \n"                       

    #y0 = [M_0,p_0,phi_0] 
    r_0 = (1./rho_0)**(.5)                        # The scaling for the radius is based on thge r^2 rho_0 = 1
    r_c = r_start / r_0                           # The starting condition needs to be rescaled 
    r_f = r_stop / r_0                            # The final radius needs o be rescaled 

    t_span = (r_c,r_f)
    t_eval = np.linspace(r_c,r_f,1000)


    M_c = 10. * 10**(-15.)                        # M = r_0 hat(M) but this is = 0 
    nu_c = 10. * 10**(-15.)                       # No scale factor for nu
    phi_c = 4. * 10**(-3.)                        # Value of scalar field at infinity 
    A4_c = (math.exp(alpha_0 * phi_c))**4.        # A(phi_c) is defined in terms of alpha and phi_c
    alpha_til_0 = alpha_0                         # Coupling constant is constant in FJBD
    psi_c = (4.*pi/3)* r_c * A4_c * alpha_til_0 * (rho_c - 3. * p_c)      

    rho_c_list.append(rho_c)

    t_eval = np.linspace(r_0,r_stop,num=10.**5.)
    

    y0 = [M_c, p_c, nu_c, phi_c, psi_c]

    #print y0
    #print y0

    #soln = solve_ivp(TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)     #I took out dense_output=True for speed
    soln = solve_ivp(STT_TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)

    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    #print soln

    M = r_0 * M 
    M = M * c**2. / G
    M = M/ M_sun
    p = p * p_0
    r = r * r_0



    radii.append(r[-1]/1000.)
    masses.append(M[-1])


    print p[0], r[-1], M[-1]



radii_gr = []
masses_gr = []


for x in pressures:


    #Initial conditions
    p_c = x                                       # The central density (In geometric units)
    #print "p_c = " + str(x) + "in g=c=1 units"
    rho_c = EOS_geometric_fromfile(p_c)           # The central presssure
    #print "rho_c = " + str(rho_c) +  "in g=c=1 units"
    rho_0 = rho_c                                 # The scale factor of density 
    p_0 = rho_0                                   # the thesis uses rho_c for both pressure and density scaling
    
    p_c = x/p_0                                   # I need to pass the code scaled values p = p_0 p hat      

    print "p_c = " + str(p_c) + " scaled \n"   
              

    

    r_c = r_start 
    r_f = r_stop 

    t_span = (r_c,r_f)
    t_eval = np.linspace(r_c,r_f,1000)


    M_c = 10. * 10**(-15.)                        
    phi_c = 10. * 10**(-15.)        
    rho_c_list.append(rho_c)

    t_eval = np.linspace(r_c,r_f,num=10.**5.)
    
    y0 = [M_c,p_c,phi_c] 

    #print y0
    #print y0

    soln = solve_ivp(TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)     #I took out dense_output=True for speed

    r = soln.t
    M = soln.y[0]
    p = soln.y[1]
    #print soln

    M = M * c**2. / G
    M = M/ M_sun
    p = p * p_0




    radii_gr.append(r[-1]/1000.)
    masses_gr.append(M[-1])


    print p[0], r[-1], M[-1]



#f= open("EOS_files/mass_radius_447.dat",'w+')
#f.write("Mass (M_sun)) \t Radius (km) \n")
"""for x in range(0,len(radii)):
    if radii[x] < 20.0:
        temp = str(masses[x]) + "\t" + str(radii[x]) + "\n"
        #f.write(temp)
        print temp"""



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

"""plt.plot(sly_ps, sly_rhos, label='sly')
plt.plot(eos_ps, eos_rhos, label='polytrope')
plt.xlabel("pressure")
plt.ylabel('density')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()"""

#plt.scatter(radii, masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())
plt.plot(radii,masses, label='STT')
plt.plot(radii_gr,masses_gr, label='GR')
plt.legend()
plt.xlabel(r'radius (m)')
plt.ylabel(r'Mass (M$_{\odot}$)')
plt.ylim(0,3)
plt.xlim(0,20)
#plt.savefig("plots/GR/mass_radius_chiraleft_447.pdf")
#plt.close()
plt.show()

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


