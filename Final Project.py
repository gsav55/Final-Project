#Turbojet Engine Problem 1
from sympy.solvers import solve
from sympy import Symbol, sqrt
import math
import matplotlib.pyplot as plt
import numpy as np

#Solves Mach Number for A/A*
#M = Symbol('M')
#solve(((1/M)*((2/(g+1)*(1+(((g-1)/2)*M**2)))**((g+1)/(2*(g-1))))-3.23,M)

# Conditions given in the problem
Pa = 18.823
Ta = 216.65
To4_max = 2300
gamma1 = 1.4
#gamma2 = 1.35
R = 287
Cp1 = (gamma1/(gamma1-1))*R/1000
#Cp2 = (gamma2/(gamma2-1))*R/1000
Ma=4.0
Fst=0.06
#hc=44000

#Efficiencies
rd=(1-0.075*(Ma-1)**(1.35))

#Setup Arrays for Graphs
Mlist = []
Ilist = []
TSFClist = []
nthlist = []
nplist = []
nolist = []
Blist = []

#Flow Conditions
Toa = Ta*(1 + ((gamma1-1)/2)*Ma**2)
print 'Toa '
print Toa
Poa = Pa*(1 + ((gamma1-1)/2)*Ma**2)**(gamma1/(gamma1-1))
print 'Poa '
print Poa
u = Ma*sqrt(gamma1*R*Ta) 
print 'u '
print u

#Inlet/Diffuser
To2=Toa
print 'To2 '
print To2
Po2=Pa*rd
print 'Po2 '
print Po2

To4 = To4_max

for Arat in np.arange(1.1,10.1,0.1):
    print '--------A/A*--------'
    print Arat
    M = Symbol('M')
    M_solve = solve(((1/M)*((2/(gamma1+1)*(1+(((gamma1-1)/2)*M**2)))**((gamma1+1)/(2*(gamma1-1))))-Arat),M)
    M4 = M_solve[1]
    print 'M4 '
    print M4
    T4 = To4/(1+((gamma1-1)/2)*M4**2)
    print 'T4 '
    print T4
    P4rat = (To4/T4)**(gamma1/(gamma1-1))
    print 'P04/P4 '
    print P4rat
    u4 = M4*sqrt(gamma1*R*T4)
    print 'u4 '
    print u4
    T2 = Symbol('T2')
    T_solve = solve(((sqrt(gamma1*R*(To2-T2))*(2/(gamma1-1))*(T4+((u4**2)/R)))/u4-(gamma1*(To2-T2)*(2/(gamma1-1))))-T2,T2)
    T2 = T_solve[1]
    print 'T2 '
    print T2

    M2 = Symbol('M2')
    M2_solve =  solve(1 + ((gamma1-1)/2)*(M2**2)-(To2/T2),M2)
    M2 = M2_solve[1]
    print 'M2'
    print M2
    P2 = Po2/((To2/T2)**(gamma1/(gamma1-1)))
    print 'P2'
    print P2
    u2 = M2*sqrt(gamma1*R*T2)
    print 'u2'
    print u2

    P4 = Symbol('P4')
    P4_solve = solve(((u2*P2/T2)/(u4*P4/T4)) - Arat,P4)
    P4 = P4_solve[0]
    print 'P4'
    print P4
    Po4 = P4*P4rat
    print 'Po4 '
    print Po4

    Po6 = Po4
    To6 = To4
    P6 = Pa

    M6 = Symbol('M6')
    M6_solve =  solve((1 + ((gamma1-1)/2)*(M6**2))**(gamma1/(gamma1-1))-(Po6/P6),M6)
    M6 = M6_solve[1]
    print 'M6'
    print M6
    T6 = To6/(1+((gamma1-1)/2)*M6**2)
    print 'T6'
    print T6
    u6 = M6*sqrt(gamma1*R*T6)
    
    

    
'''		
    I = B*(u9-u)+((1+Fb)*u7-u)
    TSFC = Fb/I
    Pav=((1+Fb)*(u7**2)/2 + B*(u9**2)/2 - (B+1)*(u**2)/2)
    Pin=Fb*hc*1000
    wp=I*u
    nth=Pav/Pin
    np=wp/Pin
    no=nth*np
    if B == B:
        print 'nth '
        print nth
        print 'Pav '
        print Pav
	
    Ilist.append([I])
    TSFClist.append([TSFC])
    nthlist.append([nth])
    nplist.append([np])
    nolist.append([no])
    Blist.append([B])
	
    
# Now to plot everything!
plt.figure(1)
plt.plot(Blist, Ilist)
plt.xlabel('Bypass Ratio, B')
plt.ylabel('Specific Thrust, I')
plt.title('I vs B')

plt.figure(2)
plt.plot(Blist, TSFClist)
plt.xlabel('Bypass Ratio, B')
plt.ylabel('TSFC')
plt.title('TSFC vs B')

plt.figure(3)
plt.plot(Blist, nthlist)
plt.xlabel('Bypass Ratio, B')
plt.ylabel('Thermal Efficiency, nth')
plt.title('Thermal Efficiency vs B')

plt.figure(4)
plt.plot(Blist, nplist)
plt.xlabel('Bypass Ratio, B')
plt.ylabel('Propulsive Efficiency, np')
plt.title('Propulsive Efficiency vs B')

plt.figure(5)
plt.plot(Blist, nolist)
plt.xlabel('Bypass Ratio, B')
plt.ylabel('Overall Efficiency, no')
plt.title('Overall Efficiency vs B')

plt.show()
'''
