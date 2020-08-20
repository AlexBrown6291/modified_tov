import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import pi, G, c, hbar, m_n, h
import math
from numpy import diff


#Constants
M_sun = 1.98892 * 10.**(30.) #SI kg
from_mev = 1.602176565 * 10.**(32)  #This is to SI, to cgs its 10^33 


#-------------------------------------------------
#  EOS FROM FILE
#-------------------------------------------------

# READ IN THE EOS 

data = np.loadtxt("EOS_files/EOS_nsat_ind1.dat",skiprows=0,delimiter='\t')
data = data.transpose()
print np.shape(data)

eos_ps = data[1]      
eos_ps = eos_ps * from_mev      
eos_Es = data[2]    
eos_Es = eos_Es * from_mev          #This is energy density
eos_rhos = eos_Es / c**2.

print "pressure = ", eos_ps
print "mass densities = ", eos_Es/c**2.


"""data = np.loadtxt("sly_compare/sly_eos_dense.dat",skiprows=1,delimiter='\t')
data = data.transpose()

eos_ps = data[1]                #pressure in cgs
eos_ps = eos_ps * 0.1           #pressure in SI
print "pressure = ", eos_ps

eos_rhos = data[2]              #mass density in cgs
eos_rhos = eos_rhos * 10**3.    #mass 
print "mass densities = ", eos_rhos"""


#-------------------------------------------------
#  SLY tables from Ingo for comparison 
#-------------------------------------------------
"""data = np.loadtxt("sly_compare/MR_sly_comp.dat",skiprows=2)
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
LAM_file = data[7]"""

#-------------------------------------------------
#  EOS Mass Radius to check curves
#-------------------------------------------------

data = np.loadtxt("EOS_files/mass_radius_1.dat",skiprows=2)
data = data.transpose()
M_file = data[0]        #Msun
R_file = data[1]        #radii

data = np.loadtxt("EOS_files/tidal_1.dat",skiprows=0)
data = data.transpose()
LAM_file = data[2]        
M_file = data[1]        #Msun
R_file = data[0]        #radii



delta = 0.

#-------------------------------------------------
#  Function definitions 
#-------------------------------------------------

def EOS_fromfile(p):
    E = np.interp(p,eos_ps,eos_rhos)

    return E


def deriv_fromfile(p):
    if p < 0:
        return 0

    p2 = p + p*delta 
    p1 = p - p*delta 
    E2 = EOS_fromfile(p2) * c**2.
    E1 = EOS_fromfile(p1) * c**2. 
    h = E2-E1
    """print "ps = ", p2, p1
    print "rhos = ", rho2, rho1
    print "diff p = ", p2-p1
    print "diff rho =", h"""
    #print h
    #print 2.*p*10**(-10.)

    deriv = 2.*p*10**(-10.)/h
    #print "deriv w/ scale= ", deriv
    #print "deriv w/o scale= ", 2.*p*10**(-10.)/h
        
    return deriv


def TOV(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]

    #print "input = ", y

    temp = r - 2. * G * M * c**(-2.)
    temp = temp * r

    temp2 = 4. * np.pi * G * r**3. * p * c**(-2.) + G * M 

    #rho = EOS(p)
    E = EOS_fromfile(p)                                   #this is energy density
    rho = E / c**2.                                       #this is mass density

    #print p
    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - (rho + p * c**(-2.)) * temp2 / temp
    dphidr = temp2 / temp 
    #print dpdr

    #print "dpdr = ", dpdr
    #print "dMdr = ", dMdr
    #print "dphidr = ",dphidr
    #print r

    if p < 0:
        dMdr = 0
        dpdr = 0
        dphidr = 0

    return [dMdr,dpdr,dphidr]

def TOV_tidal(r,y):
    M = y[0]
    p = y[1]
    phi = y[2]
    yval = y[3]

    #print "p = ", p 

    temp = r - 2. * G * M * c**(-2.)
    temp = temp * r

    temp2 = 4. * np.pi * G * r**3. * p * c**(-2.) + G * M 

    #E = EOS_fromfile(p)                                   #this is energy density
    #rho = E / c**2.                                       #this is mass density
    rho = EOS_fromfile(p)
    E = rho * c**2.
    #print p
    dMdr = 4. * np.pi * rho * r**2.
    dpdr = - (rho + p * c**(-2.)) * temp2 / temp
    dphidr = temp2 / temp 

    #print "input = ", y

    """print "Mass = ", M
    print "G = ", G
    print "r = ", r
    print "c^2= ", c**2."""

    temp1 = (2. * M * G)/ (r * c**2.)
    dpde = deriv_fromfile(p)

    F = (1. - ((4. * pi * G * r**2.)/(c**4.)) * (E - p) ) * (1. - temp1)**(-1.)


    r2Q = ((4. * pi * G * r**2.)/(c**4.)) * (5. * E + 9. * p + ((E + p)/(dpde))) * (1. - temp1)**(-1.) \
         -  6. * ( 1.- temp1)**(-1.)    \
         -  temp1**2.  * (1. + ((4. * pi * p * r**3.)/(M  * c**2.)))**2. * (1. - temp1)**(-2.)
    #print "r2Q = ", r2Q


    dydr = - (1./r) *  (yval**2. + yval * F + r2Q ) 
    #print "dydr = ", dydr
    #dydr = 100
    #print "r = ", r

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

    return [dMdr,dpdr,dphidr,dydr]

def star_boundary(r,y):
    return y[1]

#Set star boundary at pressure = 0
star_boundary.terminal = True


#-------------------------------------------------
#  Set up for integrator 
#-------------------------------------------------
#rho_0 = E_file[0]      #from Ingo's file
#p_0 = np.interp(rho_0,eos_rhos,eos_ps)
#print "Initial density = ", rho_0
#print "Initial pressure = ", p_0



phi_0 = 0. 
M_0 = 1.* 10.**-15.
r_0 = 1.* 10.**-15. #m
r_stop = 20. #km
r_stop = r_stop * 10.**(3.) #SI = m
t_span = (r_0,r_stop)
t_eval = np.linspace(r_0,r_stop,1000)
#step_size = 5. * 10**(-3.)
step_size = 1.

#p = np.logspace(25.,45.,num=10) #SI
pressures = np.logspace(33.5,35,num=10)   
pressure = 10**34.


#pressures for sly comparison
"""pressures = []
#pressures.append(p_0)
for i in range(0,len(P_c_file)):
    if i%50 == 0:
        pressures.append(P_c_file[i] * from_mev)""" 

#print len(pressures)

i = 0
radii_tidal = []
radii = []
masses = []
masses_tidal = []
lambdas = []
ks = []
steps = np.logspace(-3,1,5)
steps = steps[::-1]
print steps 

delta_p =  np.logspace(-11,-10,2)
average_error = []

for q in delta_p:
    radii_tidal = []
    radii = []
    masses = []
    masses_tidal = []
    lambdas = []
    ks = []


    for x in pressures:
        #--------------------------------------------------
        #    INTEGRATION
        #--------------------------------------------------

        # set initial conditions
        yval0 = 2.
        y0 = [M_0, x, phi_0, yval0]
        delta = q
        #y0 = [x]
        print y0
        print "rho = ", EOS_fromfile(x)

        #soln = solve_ivp(TOV_tidal,t_span,y0,method='RK45', events=star_boundary, dense_output=True)
        soln = solve_ivp(TOV_tidal,t_span,y0,method='RK45', events=star_boundary, max_step=0.1, dense_output=True)

        r = soln.t
        m = soln.y[0]
        p = soln.y[1]
        phi = soln.y[2]
        y = soln.y[3]
        #print soln

        radii_tidal.append(r[-1])
        masses_tidal.append(m[-1])

        radii.append(r[-1])
        masses.append(m[-1])

        print p[0], r[-1], m[-1]/M_sun
        #print len(M)

        #------------------------------------------------
        #Calculate the tidal love number
        #------------------------------------------------
        y_R = y[-1]
        R = r[-1]
        M = m[-1]

        beta = (G * M)/(R * c**2.)    #compactness
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
        temp = EOS_fromfile(p[x])
        #temp = temp / 10**3.
        rho_plot.append(temp)
        E_plot.append(temp*c**2.)
        temp1 = (2. * m[x] * G)/ (r[x] * c**2.)
        #print "temp1 = ", temp1
        dpde = deriv_fromfile(p[x])

        #print p[x]
        #print E_plot[x]
        #print (1. - ((4. * pi * G * r[x]**2.)/(c**4.)) * (E_plot[x] - p[x]) ) * (1. - temp1)**(-1.)
        F_plot.append((1. - ((4. * pi * G * r[x]**2.)/(c**4.)) * (E_plot[x] - p[x]) ) * (1. - temp1)**(-1.))
        r2Q_1.append(((4. * pi * r[x]**2. * G )/ c**4.) * ( 5. * E_plot[x] + 9. * p[x] + (E_plot[x]+p[x])/(dpde))*(1.-temp1)**(-1.)) # The error is here... it must come from the derivative
        r2Q_2.append(- 6. * (1. - temp1)**(-1.))
        r2Q_3.append(- temp1**2.  * (1. + ((4. * pi * r[x]**3. * p[x])/(m[x] * c**2.)))**2. * (1. - temp1)**(-2.))

        plt.plot(r/1000.,r2Q_1,color='r',label='calculated')
        plt.plot(r_file,r2q_1_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        plt.savefig("plots/GR/r2q_1_sly.pdf")
        plt.close()
        #plt.show()

        plt.plot(r/1000.,r2Q_2,color='r',label='calculated')
        plt.plot(r_file,r2q_2_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        plt.savefig("plots/GR/r2q_2_sly.pdf")
        plt.close()
        #plt.show()

        plt.plot(r/1000.,r2Q_3,color='r',label='calculated')
        plt.plot(r_file,r2q_3_file,color='b',label='from file')
        plt.xlabel(r'radius (km)')
        plt.ylabel(r'$r^{2} Q(r)_1$')
        plt.ylim(-20,50)
        plt.legend()
        plt.savefig("plots/GR/r2q_3_sly.pdf")
        plt.close()
        #plt.show()

        plt.plot(r/1000.,F_plot,color='r',label='calculated')
        plt.plot(r_file,F_file,color='b', label='from file')
        plt.xlabel("radius (km)")
        plt.ylabel("F(r)")
        plt.ylim(0,2.5)
        plt.legend()
        plt.savefig("plots/GR/f_sly.pdf")
        plt.close()
        #plt.show()

        plt.plot(r/1000.,rho_plot, label='calculated')
        plt.plot(r_file,E_file,label='from file')
        plt.yscale('log')
        plt.xlabel("radius (km)")
        plt.ylabel(r"$\rho kg/m^3$")
        plt.legend()
        plt.savefig("plots/GR/rho_sly.pdf")
        plt.close()
        #plt.ylim(0,2.5)
        plt.show()

        plt.plot(r/1000.,y, label='calculated')
        plt.plot(r_file,y_file,label='from file')
        plt.xlabel('radius (km)')
        plt.ylabel('y(r)')
        plt.legend()
        plt.savefig("plots/GR/y_sly.pdf")"""

        #plt.show()

        """plt.plot(r,M)
        plt.xlabel("radius (m)")
        plt.ylabel("pressure")()
        plt.savefig("plots/GR/mass_profile_SI_" + str(i) + ".pdf")
        plt.close()



        plt.plot(r,p)
        plt.xlim(0,r_stop)
        plt.xlabel("radius (m)")
        plt.ylabel("pressure")
        plt.savefig("plots/GR/pressure_profile_SI_" + str(i) + ".pdf")
        plt.close()
        i+=1"""

    #----------------------------------------
    # Calculate error for tidal deformability
    #----------------------------------------
    for i in range(0,len(masses_tidal)):
        #print "Radius = ", radii_tidal[x]/1000., "mass = ", masses_tidal[x]/M_sun, "\t lambda = ",  lambdas[x] #, "\t k2 = ", ks[x]
        masses_tidal[i] = masses_tidal[i]/M_sun
        radii_tidal[i] = radii_tidal[i]/1000.
        masses[i] = masses[i]/M_sun
        radii[i] = radii[i]/1000.

    avg_err = 0.

    for i in range(0,len(radii)):
        temp = np.interp(masses[i],M_file,LAM_file)
        err = np.abs((temp - lambdas[i])/temp) * 100
        avg_err = avg_err + err
        print "Radius = ", radii_tidal[i], "mass = ", masses_tidal[i], "\t lambda calculated = ",  lambdas[i], "\t lambda = ",  temp, "\t percent error = ", err

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
    plt.title(r'$P_c = 10^{34}$, $\delta P = $' + str(q) + r' step_size=0.1')
    plt.show()

print "\n"

for x in range(0,len(delta_p)):
    print "Delta = ", delta_p[x], "\t avg error = ", average_error[x]


#f= open("EOS_files/mass_radius_lambda_1436.dat",'w+')
#f.write("Mass (M_sun)) \t Radius (km) \t Dimensionless Tidal Deformability \n")

"""for x in range(0,len(masses_tidal)):
    #print "Radius = ", radii_tidal[x]/1000., "mass = ", masses_tidal[x]/M_sun, "\t lambda = ",  lambdas[x] #, "\t k2 = ", ks[x]
    masses_tidal[x] = masses_tidal[x]/M_sun
    radii_tidal[x] = radii_tidal[x]/1000.
    masses[x] = masses[x]/M_sun
    radii[x] = radii[x]/1000.
    #print masses_tidal[x], "\t", radii_tidal[x], "\n"
    #if radii_tidal[x] < 20.0:
        #temp = str(masses_tidal[x]) + "\t" + str(radii_tidal[x]) + "\t" + str(lambdas[x]) + "\n"
        #f.write(temp)"""



"""f= open("EOS_files/mass_radius_lambda_1273.dat",'w+')
f.write("Mass (M_sun)) \t Radius (km) \t Dimensionless Tidal Deformability \n")
for x in range(0,len(radii_tidal)):
    if radii[x] < 20.0:
        temp = str(masses_tidal[x]) + "\t" + str(radii_tidal[x]) + "\t" + str(lambdas[x]) + "\n"
        f.write(temp)"""



#----------------------------------------
# PLOT COMPARE CALCULATED MASS AND TIDAL TO FILES
#----------------------------------------

"""#test cgs
plt.scatter(R_file,M_file,color='r',label='from file')
plt.scatter(radii,masses, c=pressures, cmap='inferno', label='calculated')
plt.legend()
#plt.ylim(0.5,2.5)
plt.xlabel("radius (km)")
plt.ylabel(r"Mass (M$_{\odot}$)")
plt.show()

plt.scatter(R_file,LAM_file,color = 'r',label='from file')
plt.scatter(radii,lambdas, c=pressures, cmap='inferno', label='calculated')
plt.yscale('log') 
plt.legend()
plt.xlabel("radius (km)")
plt.ylabel(r"Dimensionless Tidal Deformability ($\tilde{\Lambda}$)")
plt.title(r'$P_c = 10^{34}$, $\delta P = 8 \times 10^{-11}$, step_size=0.1')
plt.savefig("plots/GR/deriv_test_6")
plt.show()"""

"""plt.scatter(steps,lambdas)
plt.xlim(10**(-6),100)
plt.xscale('log')
plt.xlabel('step size (m)')
plt.ylabel(r"$\tilde{\Lambda}$")
#plt.savefig("plots/GR/converge_test.pdf")

plt.show()"""

"""plt.scatter(delta_p,lambdas)
plt.xlim(10**(-13),10**(-6))
plt.xscale('log')
plt.xlabel(r'$\delta P$')
plt.ylabel(r"$\tilde{\Lambda}$")
plt.title(r'$P_c = 10^{34}$, $\tilde{\Lambda} \approx 4000$')
#plt.savefig("plots/GR/converge_test.pdf")
plt.show()"""

"""print radii

print "Compactness:"

for x in range(0,len(masses)):
    print "Radius = ", radii_tidal[x]/1000., "mass = ", masses_tidal[x]/M_sun, "\t lambda = ",  lambdas[x] #, "\t k2 = ", ks[x]
    masses_tidal[x] = masses_tidal[x]/M_sun
    #print "W/o Tidal calc", masses[x] * G / (c **2. * radii[x])
    masses[x] = masses[x]/M_sun"""

#radii = radii * 1000.
#plt.scatter(radii,masses,c=pressures, cmap='inferno',norm=matplotlib.colors.LogNorm())

"""plt.scatter(R_file,LAM_file,color = 'r',label='from file')
plt.scatter(radii,lambdas, c=pressures, cmap='inferno', label='calculated')
plt.yscale('log') 
plt.legend()
#plt.ylim(0.5,2.5)
plt.xlabel("radius (km)")
plt.ylabel(r"Dimensionless Tidal Deformability ($\tilde{\Lambda}$)")
#plt.show()
plt.savefig("plots/GR/sly_tidal.pdf")
plt.close()

plt.scatter(R_file,M_file,color = 'r',label='from file')
plt.scatter(radii,masses, c=pressures, cmap='inferno', label='calculated')
plt.legend()
#plt.ylim(0.5,2.5)
plt.xlabel("radius (km)")
plt.ylabel(r"Mass (M$_{\odot}$)")
#plt.show()
plt.savefig("plots/GR/sly_MR.pdf")
#plt.close()"""



