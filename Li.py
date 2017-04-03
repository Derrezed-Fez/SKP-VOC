""" 
"Robust Nonfitting Way To Determine Mass Diffusivity and Initial Concentration 
for VOCs in Building Materials with Accuracy Estimation" by Min Li.

Created by Shane Canaday
Started on 7/2/14
"""
import numpy as np

"""
Constants for this model can be seen below. Each one is capable of being
modified.
================================================================================
"""

m_length =2.00e-2                                                               #material length (m)
m_width =2.00e-2                                                                #material width(m)
m_thickness =2.54e-4                                                            #L - material layer thickness (m)
m_area = m_length * m_width                                                     #A - material area (m**2)
m_surf_area = (2 * m_length * m_width) + (2 * m_width * m_thickness) \
                + (2 * m_length * m_thickness)                                  #exposed material surface area (m**2)

c_volume = 4.45e-5                                                              #V - air chamber volume (m**3)
c_flow_rate = 1.67e-6                                                           #Q - vol air flow through chamber(m**3 s**-1)

acr = c_flow_rate/c_volume                                                      #m - air change rate (Q/V)
R = m_surf_area*m_thickness/c_volume                                            #AL/V


part_coef = 61.4397                                                             #Kv- material air partition coefficient (dimensionless)

t1 = 3600.
air_conc_t1 = 6.42404968e+02


t2 = 7200.
air_conc_t2 = 4.53346985e+02

"""
Determination of "n" using equation 7 from Min Li paper.
================================================================================
"""
n = (1/(t2-t1))*(np.log(air_conc_t1/air_conc_t2))                               #equation 7 from Min Li

#print n

"""
Determination of the root (alpha) from equation 2b from Min Li paper using
bisection method.
================================================================================
"""

def LiRoot(root):                                              
    return ((np.tan(root*m_thickness))/(root*m_thickness)) - ((acr-n)/(part_coef*R*n))#equation 2b

def bisect(LiRoot, a, b, e):                                                    #(equation, lower bound, upper bound, error tolerance)
	n = 0                                                                   #starting increment
	fa = LiRoot(a)                                                          #sets f(a) to be called with fa
	if fa == 0.0: return a                                                  #if f(LB)=0 then root=a
	fb = LiRoot(b)                                                          #sets f(b) to be called with fb
	if fb == 0.0: return b                                                  #if f(UB)=0 then root=b
	while (abs(a-b) > e):                                                   #uses error tolerance
		c = 0.5*(a+b)                                                   #defines midpoint (MP)
		fc = LiRoot(c)                                                  #sets f(c)=fc
		if fc == 0.0: return c                                          #if f(MP)=0 then root=c
		n = n + 1                                                       #increases increment by 1
		if fb*fc < 0.0:                                                 #UB*MP<0.0   
			a = c                                                   #LB becomes MP
			fa = fc	                                                #LB becomes MP
		else:
			b = c                                                   #UB becomes MP
			fb = fc                                                 #UB becomes MP
	if fa < fb:                                                             #if f(LB)<f(UB)
		return a
	else:
		return b
                                        
 
calc_root = bisect(LiRoot, 1.0e-8/m_thickness, (np.pi/2)/m_thickness, 1.0e-8)   #choosing LB, UB, and error                           

#print calc_root

"""
Determination of Diffusion Coefficient (D) using n=D*alpha**2 equation 4 from
Min Li paper.
================================================================================
"""

diff_coef = n / calc_root**2                                                    #equation 4 from Min Li

print "D = ", diff_coef

"""
Determination of "N" using equation 6 from Min Li paper.
================================================================================
"""

f1 = (acr+n)/n                                                                  #calculating f1

f2 = (acr-n)/n                                                                  #calculating f2

N = (2*R*np.tan(calc_root*m_thickness))/(calc_root*m_thickness*(f1+f2*\
    calc_root*m_thickness*((1/np.tan(calc_root*m_thickness))+\
    np.tan(calc_root*m_thickness))))                                            #equation 6 from Min LI

#print N

"""
Determination of Co using equation 8 from Min Li paper.
================================================================================
"""

Co = air_conc_t1/(N*np.exp(-1*n*t1))                                            #equation 8 Min Li

print "Co = ", Co
