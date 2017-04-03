# -*- coding: utf-8 -*-
"""
Model for the emission of volatile organic compounds

[VOC] per 30 seconds

Little, J.C.; Hodgson, A. T.; Gadgil, A. J. Modeling Emissions of Volatile 
Organic Compounds from New Carpets. Atmos. Environ. 1993, 28.

Created on Thu Jan 09 2013

Author: Stratton / Canaday
================================================================================
"""

import numpy as N
import matplotlib.pyplot as plt

"""
Constants for this model can be seen below. Each one is capable of being
modified.
================================================================================
"""
m_length = 2e-2                                                                 #length of material (m)
m_width = 2e-2                                                                  #width of material (m)
m_thickness = 2.54e-4                                                           #L - material thickness (m)
m_surf_area = (2 * m_length * m_width) + (2 * m_width * m_thickness)\
            + (2 * m_length * m_thickness)                                      #exposed material surface area (m**2)

air_volume = 4.45e-5                                                            #V - air chamber volume (m**3)
air_flow_rate = 1.67e-6                                                         #Q - vol air flow (m**3 s**-1)

diff_coef = 9.1773e-13                                                          #D - diffusion coeffient in material (m**2 s**-1)
part_coef = 61.4397                                                             #Kv- material air partition coefficient (dimensionless)
start_conc = 1.4464e8                                                           #Co - initial material [VOC] (ug m**-3)

h = (air_flow_rate / m_surf_area) / (diff_coef * part_coef)
k = (air_volume / m_surf_area) / part_coef

time = N.zeros((481,1), dtype='f')                                              #t - time (s)
delta_time = 30                                                                 #delta t - time increment (s)
lin_dist = m_thickness                                                          #x - was multiplied by a loading ratio of 1
dx = 0
tmp_root = 0

"""
Definitions of root equation from paper (littleroot) and oddNum used to 
generate a list of odd numbers used for index numbers (see Excel spreadsheet).
================================================================================
"""
root = N.zeros((len(time),1), dtype='f')                                        #creates an array of zeroes 1 column by xxx rows

def littleRoot(root, m_thickness, h, k):                                        #defining littleroot
	littleRoot = N.tan(root * m_thickness) - (h / root) + (k * root)        #equation 10 from Little et al. 
	return littleRoot
 
def oddNum(start):                                                              ##start: value to begin returning odd numbers from
    while True:
        if start % 2 == 0:                                                      #start value / d has a remainder = 0
            start += 1
        yield start
        start += 1

someOddNums = oddNum(1)                                                         #generates odd numbers using above function

bounds = N.zeros((len(time), 3), dtype='f')                                     #initialize an array of lenth(time)row by 3 col wide for index [0], lower [1] and
                                                                                #upper bounds [2]
calc_root = N.zeros((len(time),1), dtype='f')

air_conc = N.zeros((len(time),1), dtype='f')

time = N.arange(0, delta_time * len(time), delta_time, dtype = 'f')             #N.arange([start], stop[, step], dtype=None)
                                                                                #returns evenly spaced values within the interval step is spacing between values

"""
Uses Bisection method to determine roots (f(x)=0). Creates upper (a) and lower 
(b) bounds and then takes midpoint (m=b-a/2).The function is evaluated and the 
signs of f(a), f(b), f(m) are compared.
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
    for j in xrange(len(time)):                                                 #sum portion of little equation
        series_sum += (N.exp(-1 * diff_coef * calc_root[j]**2 * time[i]) *      #equation 7 from Little et al.
        (h - k * calc_root[j]**2) * N.cos(calc_root[j] * lin_dist)) / \
        ((m_thickness * (h - k * calc_root[j]**2)**2 + calc_root[j]**2 * 
        (m_thickness + k) + h) * N.cos(calc_root[j] * m_thickness))
        
    air_conc[i] = (2 * start_conc * series_sum) / part_coef                     #equation 5 from Little et al.

"""
Experimental Data points to be plotted on same graph
================================================================================
"""
#exp_time = N.array([0,3600,7200,10800,14400,18000,21600,25200,28800,32400,36000,39600,43200,46800,50400,54000,57600,61200,64800,68400,72000,75600,79200,82800,87480,91080,94680,98280,101880,105480,109080,112680,116280,119880,123480,127080,130680,134280,137880,141480,145080,148680,152280,155880,159480,163080,166680,170280,173880,177480,181080,184680,188280,191880,195480,199080,202680,206280,209880,213480,217080,219720,222480,225360,228180,231480,235080,238680,242280,245880,251760,255360,258960,262560,266160,269760,273360,276960,280560,284160,287760,291360,294960,298560,302160,305760,309360,312960,316560,320160,323760,327360,330960,334560,338160])                                                     # Create a step function (will demonstrate monotonicity)
#exp_conc = N.array([0,119.3304286,166.8987075,127.3791435,114.5518041,89.67933297,79.24893821,69.68047201,60.30824896,56.0943794,51.81676047,38.91206524,36.24089355,35.91131135,33.04574829,34.24561649,34.07034861,28.9155855,29.33526387,28.99588238,25.49463646,25.44901334,24.67863796,24.95762622,21.47983346,22.36881623,21.60750922,21.10756538,21.11481666,18.90938089,18.1525141,19.65910815,17.09027817,17.0965776,15.77917156,14.68400953,14.88755069,14.68553593,15.9505889,14.93721037,14.96502871,13.06755316,14.65597137,14.94226277,13.51866212,13.03338354,13.28472431,12.12041203,12.32173917,11.27170754,11.73942396,10.82070051,10.39801253,10.75020939,10.22978854,10.44517916,10.90296822,9.225371582,9.096878121,9.484838614,9.02282305,9.839133936,9.637473537,8.252501731,8.669481861,8.81452524,8.69173416,8.399208818,8.291565076,7.380819218,7.340743937,7.544517317,7.142052607,6.443223373,6.366104866,6.902176708,6.995562006,6.308661821,5.912827382,5.954208353,6.195959124,5.998746578,5.705085047,6.110455583,5.888671752,5.671292433,5.558686038,5.300870919,5.24881991,4.548396346,4.688283491,4.781920546,4.880534247,4.591161187,4.549769979])

"""
Creation of Plot
================================================================================
"""
plt.plot(time / 60, air_conc, 'b')                                              #model prediction curve
#plt.plot(exp_time / 60, exp_conc, 'ro')                                        #experimental data curve
plt.title('[VOC] vs Time')
plt.xlabel('$Time\ (min)$')
plt.ylabel('$Concentration\ (\mu g/m^3)$')
legends = ["Model", "Data"]
plt.legend(legends, loc="upper right")
plt.show()

#print time
#print air_conc
