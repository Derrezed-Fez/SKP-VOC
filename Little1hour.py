# -*- coding: utf-8 -*-
"""
Model for the emission of volatile organic compounds

[VOC] per hour

Little, J.C.; Hodgson, A. T.; Gadgil, A. J. Modeling Emissions of Volatile 
Organic Compounds from New Carpets. Atmos. Environ. 1993, 28.

Created on Thu Jan 09 2013

Author: Stratton / Canaday
================================================================================
"""

import numpy as N
import matplotlib.pyplot as plt

"""
Parameters for this model can be seen below. Each one is capable of being
modified.
================================================================================
"""
m_length =0.10                                                                  #length of material (m)
m_width =0.10                                                                   #width of material (m)
m_thickness =2.5e-4/2                                                           #L - material layer thickness (m)
m_surf_area = (2 * m_length * m_width) #+ (2 * m_width * m_thickness)\
              #+ (2 * m_length * m_thickness)                                   # exposed material surface area (m**2)

air_volume = 0.053                                                              #V - air chamber volume (m**3)
air_flow_rate = 1.472e-5                                                        #Q - vol air flow through chamber(m**3 s**-1)

diff_coef = 1.736e-13                                                           #D - diffusion coeffient in material (m**2 s**-1)
part_coef = 230                                                                 #Kv- material air partition coefficient (dimensionless)
start_conc = 1.566e8                                                            #Co - initial material [VOC] (ug m**-3)

h = (air_flow_rate / m_surf_area) / (diff_coef * part_coef)
k = (air_volume / m_surf_area) / part_coef

time = N.zeros((148,1), dtype='f')                                              #t - time (s)
delta_time = 3600                                                               #delta t - time increment (s)
lin_dist = m_thickness                                                          #x - was multiplied by a loading ratio of 1
dx = 0
tmp_root = 0

"""
Definitions of root equation from paper (littleroot) and oddNum used to 
generate a list of odd numbers used for index numbers (see Excel spreadsheet).
================================================================================
"""
root = N.zeros((len(time),1), dtype='f')                                        #creates an array of zeroes 1 column by length(time) rows

def littleRoot(root, m_thickness, h, k):                                        #defining littleroot 
	littleRoot = N.tan(root * m_thickness) - (h / root) + (k * root)        #equation 10 from Little et al.
	return littleRoot
 
def oddNum(start):                                                              #start: value to begin returning odd numbers from
    while True:
        if start % 2 == 0:                                                      #start value / d has a remainder = 0
            start += 1
        yield start
        start += 1

someOddNums = oddNum(1)                                                         #generates odd numbers using above function

bounds = N.zeros((len(time), 3), dtype='f')                                     #initialize an array of length(time) rows by 3 columns wide for index [0], lower [1] and
                                                                                #upper bounds [2]
calc_root = N.zeros((len(time),1), dtype='f')

air_conc = N.zeros((len(time),1), dtype='f')

time = N.arange(0, delta_time * len(time), delta_time, dtype = 'f')             #N.arange([start], stop[, step], dtype=None)
                                                                                #returns evenly spaced values within the interval step is spacing between values
time.shape = (len(time),1)                                                      #sets time.shape to (shapebounds [0], 1) or (101,1)

"""
Uses Bisection method to determine roots (f(x)=0). Creates upper (a) and lower 
(b) bounds and then takes midpoint (m=b-a/2).The function is evaluated and the 
signs of f(a), f(b), f(m) are compared
================================================================================
"""
for i in xrange(len(time)):
    bounds[i][0] = someOddNums.next()                                           #integers for index ([0]) loops through odd numbers making column 1
    bounds[i][2] = bounds[i][0] * (N.pi / (2 * m_thickness))                    #calculates upper bounds and puts in column 3
    bounds[i][1] = bounds[i - 1][2]                                             #calculates lower bounds and puts in column 2
for i in xrange(len(time)):
    lb = bounds[i][1]                                                           #sets lb (lower bounds) to column 2 of array
    ub = bounds[i][2]                                                           #sets ub (upper bounds) to column 3 of array 
    root[i] = lb
    dx = ub - lb
    for j in xrange(500):
        dx = dx * 0.5                                                           #midpoint calculation
        calc_root[i] = root[i] + dx                                             #calc_root[i] = QN of excel spreadsheet
        temp_root = littleRoot(calc_root[i], m_thickness, h, k)                 #temp_root = FQN in excel spreadsheet
        if temp_root < 0:
            root[i] = calc_root[i]
        if abs(dx) < 1e-100 or temp_root == 0: 
            break

"""
Calculate the concentration of VOC in the chamber air and fill in the 
air_conc array.
================================================================================
"""
for i in xrange(len(time)):
    series_sum = 0                                                              #set intial sum of little equation to zero
    for j in xrange(len(time)):                                                 #sum portion of Little equation
        series_sum += (N.exp(-1 * diff_coef * calc_root[j]**2 * time[i]) *      #equation 7 from Little et al.
        (h - k * calc_root[j]**2) * N.cos(calc_root[j] * lin_dist)) / \
        ((m_thickness * (h - k * calc_root[j]**2)**2 + calc_root[j]**2 * 
        (m_thickness + k) + h) * N.cos(calc_root[j] * m_thickness))
        
    air_conc[i] = (2 * start_conc * series_sum) / part_coef                     #equation 5 from Little et al.

"""
Experimental Data points to be plotted on same graph
================================================================================
"""
exp_time = N.array([396.0,1188.0,1800.0,2520.0,3600.0,7200.0,14832.0,21600.0,28800.0,36000.0,86400.0,100800.0,115200.0,172800.0,389088.0,441360.0,525996.0])                                                     
exp_conc = N.array([527.19,696.7,803.92,825.94,929.58,886.83,740.58,564.39,411.42,307.23,60.42,40.32,32.58,13.96,5.97,4.47,4.93])

"""
Creation of Plot
================================================================================
"""
plt.plot(time / 3600, air_conc, 'b')                                            #model prediction curve
plt.plot(exp_time / 3600, exp_conc, 'ro')                                       #experimental data curve
plt.title('[VOC] vs Time')
plt.xlabel('$Time\ (hr)$')
plt.ylabel('$Concentration\ (\mu g/m^3)$')
legends = ["Model", "Data"]
plt.legend(legends, loc="upper right")
plt.show()

print air_conc