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



# READ IN THE EOS 


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


data = np.loadtxt("sly_compare/sly_eos_dense.dat",skiprows=1,delimiter='\t')
data = data.transpose()
eos_ps = data[1]                #pressure in cgs
eos_ps = eos_ps * 0.1           #pressure in SI
#print "pressure = ", eos_ps
eos_pg = eos_ps * G * c**(-4.)              #conversion to dimensionless quantities


eos_rhos = data[2]              #mass density in cgs
eos_rhos = eos_rhos * 10**3.    #mass 
#print "mass densities = ", eos_rhos
eos_rhog = eos_rhos * G * c**(-2.)          #Use G * c**-2 because this is mass density

#-------------------------------------------------
#  SLY tables from Ingo for comparison 
#-------------------------------------------------
data = np.loadtxt("sly_compare/MR_sly_comp.dat",skiprows=2)
data = data.transpose()

r_file = data[0]                # km
E_file = data[1] * 10**3.       #kg/m^3
y_file = data[2]
F_file = data[3]
r2q_file = data[4]
Cs_file = data[5]
r2q_1_file = data[6]
r2q_2_file = data[7]
r2q_3_file = data[8]

data = np.loadtxt("sly_compare/MR_sly.dat",skiprows=2)
data = data.transpose()

R_file = data[0]        #km
M_file = data[1]        #M_sun
P_c_file = data[2]
rh0_c_file = data[3]   
beta_file = data[5] 
k2_file = data[6]
LAM_file = data[7]




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


def deriv_fromfile(p):
    if p < 0:
        return 0
    P = p * p_c
    p2 = P + P*10**(-10.)
    p1 = P - P*10**(-10.)
    E2 = EOS_geometric_fromfile(p2) 
    E1 = EOS_geometric_fromfile(p1) 
    h = E2-E1
    """print "ps = ", p2, p1
    print "rhos = ", rho2, rho1
    print "diff p = ", p2-p1
    print "diff rho =", h"""
    #print h
    #print 2.*p*10**(-10.)

    deriv = 2.*P*10**(-10.)/h
    #print "deriv w/ scale= ", deriv
    #print "deriv w/o scale= ", 2.*p*10**(-10.)/h
        
    return deriv


def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    rho = EOS_scaled_fromfile(p)
    
    temp1 = r * (r - 2. * M )

    temp2 = M + 4 * np.pi * r**3.  * p * p_c 

    dMdr = 4. * np.pi * rho * rho_c * r**2.
    dpdr = - 1. * (p * p_c + rho * rho_c) * ( temp2 / temp1 ) * (1./p_c)
    dphidr = temp2 / temp1 



    return [dMdr,dpdr,dphidr]

def TOV_tidal(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]
    yval = y[3]

    rho = EOS_scaled_fromfile(p)
    E = rho
    
    temp1 = r * (r - 2. * M )
    temp2 = M + 4 * np.pi * r**3.  * p * p_c 

    dMdr = 4. * np.pi * rho * rho_c * r**2.
    dpdr = - 1. * (p * p_c + rho * rho_c) * ( temp2 / temp1 ) * (1./p_c)
    dphidr = temp2 / temp1 


    temp1 = (2. * M)/ (r)
    dpde = deriv_fromfile(p)

    F = (1. - ((4. * pi * r**2.) * (E * E_c  - p * p_c) )) * (1. - temp1)**(-1.)


    r2Q = (4. * pi * r**2.) * (5. * E * E_c + 9. * p * p_c + ((E * E_c + p * p_c)/(dpde))) * (1. - temp1)**(-1.) \
         -  6. * ( 1.- temp1)**(-1.)    \
         -  temp1**2.  * (1. + ((4. * pi * p * p_c * r**3.)/(M)))**2. * (1. - temp1)**(-2.)
    #print "r2Q = ", r2Q


    dydr = - (1./r) *  (yval**2. + yval * F + r2Q ) 

    if M == 0:
        dydr = 0


    #dydr = yval

    #print "f = ", F
    #print "r2q = ", r2Q
    #print "dydr = ", dydr

    if p < 0:
        dMdr = 0
        dpdr = 0
        dphidr = 0
        dydr = 0

    """print "dpdr = ", dpdr
    print "dMdr = ", dMdr
    print "dphidr = ",dphidr
    print "dydr =", dydr"""
    #dydr = 0
    return [dMdr,dpdr,dphidr,dydr]




def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True

 
M_0 = 10. * 10**(-15.)
r_0 = 0.0001 #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)

p_low = 1. * from_mev                           #Ingo says to start with ~ 1MeV cm^-3
#print p_low
sly_density = E_file[0] * G * c**(-2.)

pressures = np.logspace(33.5,35,num=10)           #p_max in the file from Ingo is 2.4 x 10^35, so max needs to be less than that 
#pressures = []
#pressures.append(np.interp(sly_density,eos_rhog,eos_pg))
pressures = pressures * G * c**(-4.)
print pressures
#print pressures



radii = []
masses = []
lambdas = []
ks = []

nums = np.linspace(10.**5.,10.**8.,num=100)
#step_size = np.linspace(1.,10**(-5.))
step_size = 1.
#print step_size
rho_c_list = []
steps = np.logspace(-3,1,5)
steps = steps[::-1]
print steps 


average_error = []

for q in steps:
    radii = []
    masses = []
    lambdas = []
    ks = []

    for x in pressures:

        #Initial conditions
        yval = 2.
        p_c = x
        rho_c = EOS_geometric_fromfile(p_c)
        E_c = rho_c
        rho_c_list.append(rho_c)

        t_eval = np.linspace(r_0,r_stop,num=10.**5.)
        p_0 = x/p_c
        phi_0 = 0.
        y0 = [M_0,p_0,phi_0,yval]
        print y0
        #print y0
        #y0 = [x]
        #print y0

        #soln = solve_ivp(TOV, t_span, y0,  method='RK45', events=star_boundary, max_step=step_size)     #I took out dense_output=True for speed
        soln = solve_ivp(TOV_tidal,t_span,y0,method='RK45', events=star_boundary, max_step=q, dense_output=True)
        r = soln.t
        m = soln.y[0]
        p = soln.y[1]
        y = soln.y[3]


        #Calculate the tidal love number
        #------------------------------------------------
        y_R = y[-1]
        R = r[-1]
        M = m[-1]

        beta = (M)/(R)    #compactness
        tmp = 1. - 2* beta

        k2 = ((8. * beta**5.)/5.) * (1. - 2. * beta)**2. * (2. - y_R + (y_R - 1.) * 2. * beta)  \
            * ( 2. * beta * (6. - 3. * y_R + 3. * beta * (5. * y_R - 8.)) \
            + 4. * beta**3. * (13. - 11. * y_R + beta * (3.* y_R - 2.) + 2. * beta**2. * (1. + y_R)) \
            + 3. * (1. - 2. * beta)**2. * (2. - y_R + 2. * beta * (y_R - 1.)) * np.log(1. - 2. * beta)) **(-1.)

        #print "k_2 (one) ", k2



        l_dim = 2. * k2 * beta**(-5.) / 3
        ks.append(k2)
        lambdas.append(l_dim)

        
        #------------------------------------------------
        # Diagnostic  Plots
        #------------------------------------------------
        """rho_plot = []
        E_plot = []
        F_plot = []
        r2Q_1 = []
        r2Q_2 = []
        r2Q_3 = []


        for x in range(0,len(p)):
            temp = EOS_scaled_fromfile(p[x])

            rho_plot.append(temp * c**2. / G)
            E_plot.append(temp)
            temp_new = (2. * m[x])/ (r[x])
            #print "temp1 = ", temp1

            dpde = deriv_fromfile(p[x])
            F_plot.append((1. - ((4. * pi * r[x]**2.) * (temp * E_c  - p[x] * p_c) )) * (1. - temp_new)**(-1.))
            r2Q_1.append((4. * pi * r[x]**2.) * (5. * temp * E_c + 9. * p[x] * p_c + ((temp * E_c + p[x] * p_c)/(dpde))) * (1. - temp_new)**(-1.))
            r2Q_2.append(-  6. * ( 1.- temp_new)**(-1.) )   
            r2Q_3.append(-  temp_new**2.  * (1. + ((4. * pi * p[x] * p_c * r[x]**3.)/(m[x])))**2. * (1. - temp_new)**(-2.))



            rho_plot[x] = rho_plot[x] * c**2. / G

        plt.plot(r/1000.,r2Q_1,color='r',label='calculated')
        plt.plot(r_file,r2q_1_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        #plt.savefig("plots/GR/r2q_1_sly.pdf")
        #plt.close()
        plt.show()

        plt.plot(r/1000.,r2Q_2,color='r',label='calculated')
        plt.plot(r_file,r2q_2_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        #plt.savefig("plots/GR/r2q_2_sly.pdf")
        #plt.close()
        plt.show()

        plt.plot(r/1000.,r2Q_3,color='r',label='calculated')
        plt.plot(r_file,r2q_3_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        #plt.savefig("plots/GR/r2q_3_sly.pdf")
        #plt.close()
        plt.show()

        plt.plot(r/1000.,F_plot,color='r',label='calculated')
        plt.plot(r_file,F_file,color='b', label='from file')
        plt.xlabel("radius (km)")
        plt.ylabel("F(r)")
        plt.ylim(0,2.5)
        plt.legend()
        #plt.savefig("plots/GR/f_sly.pdf")
        #plt.close()
        plt.show()
        
        plt.plot(r/1000.,rho_plot, label='calculated')
        plt.plot(r_file,E_file,label='from file')
        plt.yscale('log')
        plt.xlabel("radius (km)")
        plt.ylabel(r"$\rho kg/m^3$")
        plt.legend()
        #plt.savefig("plots/GR/rho_sly.pdf")
        #plt.close()
        #plt.ylim(0,2.5)
        plt.show()

        plt.plot(r/1000.,y, label='calculated')
        plt.plot(r_file,y_file,label='from file')
        plt.xlabel('radius (km)')
        plt.ylabel('y(r)')
        plt.legend()
        plt.show()
        #plt.savefig("plots/GR/y_sly.pdf")"""
        #print soln

        M = M * c**2. / G
        M = M/ M_sun
        #p = p * p_c



        radii.append(r[-1]/1000.)
        masses.append(M)


    #----------------------------------------
    # Calculate error for tidal deformability
    #----------------------------------------
    """for i in range(0,len(masses)):
        masses[i] = masses[i]/M_sun
        radii[i] = radii[i]/1000."""

    avg_err = 0.

    for i in range(0,len(radii)):
        temp = np.interp(masses[i],M_file,LAM_file)
        err = np.abs((temp - lambdas[i])/temp) * 100
        avg_err = avg_err + err
        print "Radius = ", radii[i], "mass = ", masses[i], "\t lambda calculated = ",  lambdas[i], "\t lambda = ",  temp, "\t percent error = ", err

    avg_err = avg_err/len(radii)

    print "average error = ", avg_err
    average_error.append(avg_err)

    """plt.scatter(R_file,M_file,color='r',label='from file')
    plt.scatter(radii,masses, c=pressures, cmap='inferno', label='calculated')
    plt.legend()
    #plt.ylim(0.5,2.5)
    plt.xlabel("radius (km)")
    plt.ylabel(r"Mass (M$_{\odot}$)")
    plt.show()"""

    plt.scatter(R_file,LAM_file,color = 'r',label='from file')
    plt.scatter(radii,lambdas, c=pressures, cmap='inferno', label='calculated')
    plt.yscale('log') 
    plt.legend()
    plt.xlabel("radius (km)")
    plt.ylabel(r"Dimensionless Tidal Deformability ($\tilde{\Lambda}$)")
    plt.title(r' step_size=' + str(q))
    #plt.savefig('/work/stephanie.brown/WWW/alternate_theories/tidal/stepsize_test_'+str(q)+'.png')
    plt.show()
    plt.close()

print "\n"

for x in range(0,len(delta_p)):
    print "stepsize = ", steps[x], "\t avg error = ", average_error[x]


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


"""plt.scatter(R_file,M_file,color='k')
#plt.scatter(radii,masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())
plt.scatter(radii,masses,color='r')
plt.legend()
plt.xlabel(r'radius (m)')
plt.ylabel(r'Mass (M$_{\odot}$)')
#plt.xlim(0,20)
plt.ylim(0,3)
#plt.savefig("plots/GR/mass_radius_chiraleft_447.pdf")
#plt.close()
plt.show()

plt.scatter(R_file,LAM_file,color = 'k',label='from file')
plt.scatter(radii,lambdas, c=pressures, cmap='inferno', label='calculated')
plt.yscale('log') 
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel(r"Dimensionless Tidal Deformability ($\tilde{\Lambda}$)")
plt.title(r'$step_size=0.1')
#plt.savefig("plots/GR/deriv_test_6")
plt.show()"""

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


