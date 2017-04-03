# -*- coding: utf-8 -*-
"""
Model for the Emission of Volatile Organic Compounds (VOCs) in a Static
(Air-Tight) Environment.

Xiong, J.; Liu, C.; Zhang, Y. A General Analytical Model for Formaldehyde 
and VOC Emission/Sorption in Single-Layer Building Materials and Its Application
in Determining the Characteristic Parameters. Atmos. Environ. 2012, 47, 288–294.

Created on Thu July 24 2014

Author: Canaday
================================================================================
"""

import numpy as N
import matplotlib.pyplot as plt

"""
Parameters for this model can be seen below. Each one is capable of being
modified.
================================================================================
"""

m_length = 0.2                                                                  #material length (m)
m_width = 0.2                                                                   #material width (m)
m_thickness = 0.0084                                                            #material thickness (m)
m_surf_area = (2 * m_length * m_width) + (2 * m_width * m_thickness) \
                + (2 * m_length * m_thickness)                                  #exposed material surface area (m**2)

air_volume = 0.03                                                               #air volume of chamber (m**3)
air_velocity = 0.138                                                            #velocity of air over material (m/s)
kin_visc_air = 6.51e-7                                                          #kinematic viscosity of air (m**2/s)
air_flow_rate = 0.0                                                             #vol air flow through chamber(m**3 s**-1)

start_conc = 5.02e8                                                             #initial [VOC] in material (ug/m**3)
part_coef = 920                                                                 #air-material partition coefficient
diff_coef = 1.04e-10                                                            #diffusion coefficient for material (m**2/s)

time = N.zeros((95,1), dtype='f')                                               #length of experiment (number of delta times)                               
delta_time = 3600                                                               #delta time (sec)

"""
Convective Mass Transfer Coefficient Calculation (hm)
================================================================================
"""
Reynolds = (air_velocity*m_thickness)/kin_visc_air

Schmidt = kin_visc_air/diff_coef

Sherwood = 0.664*(Schmidt**(1/3))*(Reynolds**(1/2))

hm = (Sherwood*diff_coef)/m_thickness

"""
Alpha, Beta, and Bim Calculation
================================================================================
"""
alpha = (air_flow_rate*(m_thickness**2))/(diff_coef*air_volume)

beta = (m_surf_area*m_thickness)/air_volume

Bim = (hm*m_thickness)/diff_coef

"""
Determination of the root (qn) using bisection method.
================================================================================
"""
root = N.zeros(len(time), dtype='f')

lin_dist = m_thickness                                                          #x - was multiplied by a loading ratio of 1
dx = 0
tmp_root = 0

def oddNum(start):                                                              #start: value to begin returning odd numbers from
    while True:
        if start % 2 == 0:                                                      #start value / d has a remainder = 0
            start += 1
        yield start
        start += 1

someOddNums = oddNum(1)                                                         #generates odd numbers using above function

bounds = N.zeros((len(time), 3), dtype='f')                                     #initialize an array of 101 rows by 3 columns wide for index [0], lower [1] and
                                                                                #upper bounds [2]
calc_root = N.zeros((len(time),1), dtype='f')

time = N.arange(0, delta_time * len(time), delta_time, dtype = 'f')             #N.arange([start], stop[, step], dtype=None)
                                                                                #returns evenly spaced values within the interval step is spacing between values
time.shape = (len(time),1)  

def XiongRoot(root, alpha, part_coef, beta, Bim):                                              
    return ((alpha-(root**2))/((part_coef*beta)+(alpha-(root**2))*(part_coef \
           *(Bim**(-1))))) - (root*N.tan(root))

for i in xrange(len(time)):
    bounds[i][0] = someOddNums.next()                                           #integers for index ([0]) loops through odd numbers making column 1
    bounds[i][2] = bounds[i][0] * (N.pi / (2 * m_thickness))                    #calculates upper bounds and puts in column 3
    bounds[i][1] = bounds[i - 1][2]                                             #calculates lower bounds and puts in column 2 using previous upperbound
for i in xrange(len(time)):
    lb = bounds[i][1]                                                           #sets lb (lower bounds) to column 2 of array
    ub = bounds[i][2]                                                           #sets ub (upper bounds) to column 3 of array 
    root[i] = lb
    dx = ub - lb
    for j in xrange(500):
        dx = dx * 0.5                                                           #midpoint calculation
        calc_root[i] = root[i] + dx                                             #calc_root[i] = QN of excel spreadsheet
        temp_root = XiongRoot(calc_root[i], alpha, part_coef, beta, Bim)                 #temp_root = FQN in excel spreadsheet
        if temp_root < 0:
            root[i] = calc_root[i]
        if abs(dx) < 1e-100 or temp_root == 0: 
            break

#print calc_root
            
"""
Calculation of Gn
================================================================================
"""
calc_Gn = N.zeros((len(time),1), dtype ='f')

def Gn(calc_root):
    C1 = part_coef * beta
    C2 = alpha - (calc_root**2)
    C3 = part_coef * (Bim**(-1))
    C123 = C1+(C2 * C3)+2
    D1 = calc_root**2
    D2 = N.cos(calc_root)
    D  = C123 * D1 * D2
    E1 = alpha - (3 * D1)
    E  = C1 + (E1 * C3) + C2
    F1 = N.sin(calc_root)
    F  = E * calc_root * F1
    return D + F
    """
    return (((part_coef*beta)+((alpha-(calc_root**2))*(part_coef*(Bim**(-1))))+2)\
    *((calc_root**2)*N.cos(calc_root)))+(((part_coef*beta)+(alpha-(3* \
    (calc_root**2)))*(part_coef*(Bim**(-1))))+alpha-(calc_root**2))* \
    (calc_root*N.sin(calc_root))
    """
for i in xrange(len(time)):
    calc_Gn[i] = Gn(calc_root[i])

#print calc_Gn

"""
Calculate the concentration of VOC in the chamber air and fill in the 
air_conc array.
================================================================================
"""
air_conc = N.zeros((len(time),1), dtype='f')
calc_sum = N.zeros((len(time),1), dtype='f')
def series_sum(calc_root,calc_Gn,time):
    A1 = beta
    A2 = calc_root
    A3 = N.sin(calc_root)
    A123 = A1*A2*A3
    A4 = 0.0
    A5 = N.cos(calc_root)
    A45 = A4*A5
    A = (A123 - A45)/calc_Gn
    B1 = -1.0*diff_coef
    B2 = m_thickness**(-2)
    B3 = calc_root**2
    B4 = time 
    B1234 = B1*B2*B3*B4
    B = N.exp(B1234)
    return (sum( A * B ) for i in xrange(len(time)))

"""
for i in xrange(len(time)):
    series_sum = 0                                                              
    for j in xrange(len(time)):                                           
        
        A1 = beta
        A2 = calc_root[j]
        A3 = N.sin(calc_root[j])
        A123 = A1*A2*A3
        A4 = 0.0
        A5 = N.cos(calc_root[j])
        A45 = A4*A5
        A = (A123 - A45)/calc_Gn[j]
        B1 = -1.0*diff_coef
        B2 = m_thickness**(-2)
        B3 = calc_root[j]**2
        B4 = time[i] 
        B1234 = B1*B2*B3*B4
        B = N.exp(B1234)
        series_sum += ( A * B )

        series_sum += ((beta*calc_root[j]*N.sin(calc_root[j]-(0.0* \
        N.cos(calc_root[j]))))/calc_Gn[j])*(N.exp(-1*diff_coef*(m_thickness**2)\
        *(calc_root[j]**2)*time[i]))
    """
for i in xrange(len(time),1):  
    air_conc[i] = ((start_conc*beta)/((part_coef*beta)+1))+(2 * start_conc * series_sum(calc_root[i],calc_Gn[i],time[i]))  
    
print air_conc