"""
Determination of Co and D using:
Wei Ye, John C. Little ,Doyun Won c, Xu Zhang. Screening-Level Estimates of 
Indoor Exposure to Volatile Organic Compounds Emitted From Building Materials.
Build. Eviron. 2014 75, 58-66.

Author: Canaday
2014/07/07
================================================================================
"""

import pylab as P
from scipy.integrate import simps
"""
Parameters
================================================================================
"""
m_length = 0.10                                                                 # material length (m)
m_width = 0.10                                                                  # material width (m)
m_thickness = 2.5e-4/2                                                           # material thickness (m)
m_surf_area = (2 * m_length * m_width) #+ (2 * m_width * m_thickness) \
            #+ (2 * m_length * m_thickness)                                      # material surface area (m2)                     

air_flow_rate = 1.472e-5                                                         # air flow rate (m3/s)

peak_emis_time = 3600.0                                                         # peak emission experimental time point (s)
last_exp_time = 525996.0                                                        # last experimental time point (s)

exp_time = P.array([0.0,396.0,1188.0,1800.0,2520.0,3600.0,7200.0,14832.0,21600.0,28800.0,36000.0,86400.0,100800.0,115200.0,172800.0,389088.0,441360.0,525996.0,720000.0])                                              
exp_conc = P.array([0.0,527.19,696.7,803.92,825.94,929.58,886.83,740.58,564.39,411.42,307.23,60.42,40.32,32.58,13.96,5.97,4.47,4.93,0.0])
"""
Determination of Diffusion Coefficient (D) using Ye method:
================================================================================
"""

D_lower = (0.2  * (m_thickness)**2) / last_exp_time                             # diff. coef. lower bound
D_upper = (0.01 * (m_thickness)**2) / peak_emis_time                            # diff. coef. upper bound

print D_lower, "<= D <=", D_upper

"""
Determination of Completed Time Upper (tu) and Lower (tl) bounds:
================================================================================
"""

comp_time_lower = (2.0 * (m_thickness)**2) / D_upper
comp_time_upper = (2.0 * (m_thickness)**2) / D_lower

print "tlower = ",comp_time_lower
print "tupper = ",comp_time_upper

"""
Piecewise cubic Hermite inter-polating polynomial code was obtained from: 
http://matplotlib.1069221.n5.nabble.com/quot-Piecewise-Cubic-Hermite-
Interpolating-Polynomial-quot-in-python-td24157.html

Piecewise cubic Hermite interpolation (monotonic...) in Python

Author: Michalski
================================================================================
"""
def pchip(x, y, xnew):

    m    = pchip_init(x, y)                                                     # Compute the slopes used by the piecewise cubic Hermite interpolator
    ynew = pchip_eval(x, y, xnew)                                               # Use these slopes (along with the Hermite basis function) to interpolate          

    return ynew

#===============================================================================
def x_is_okay(x,xvec):

    n = len(x)                                                                  # Make sure "x" and "xvec" satisfy the conditions for
    m = len(xvec)                                                               # running the pchip interpolator

    xx = x.copy()                                                               # Make sure "x" is in sorted order (brute force, but works...)
    xx.sort()
    total_matches = (xx == x).sum()
    if total_matches != n:
        print "*" * 50
        print "x_is_okay()"
        print "x values weren't in sorted order --- aborting"
        return False

    delta = x[1:] - x[:-1]                                                      # Make sure 'x' doesn't have any repeated values
    if (delta == 0.0).any():
        print "*" * 50
        print "x_is_okay()"
        print "x values weren't monotonic--- aborting"
        return False

    check = xvec > x[-1]                                                        # Check for in-range xvec values (beyond upper edge)
    if check.any():
        print "*" * 50
        print "x_is_okay()"
        print "Certain 'xvec' values are beyond the upper end of 'x'"
        print "x_max = ", x[-1]
        indices = P.compress(check, range(m))
        print "out-of-range xvec's = ", xvec[indices]
        print "out-of-range xvec indices = ", indices
        return False

    check = xvec< x[0]                                                          # Second - check for in-range xvec values (beyond lower edge)
    if check.any():
        print "*" * 50
        print "x_is_okay()"
        print "Certain 'xvec' values are beyond the lower end of 'x'"
        print "x_min = ", x[0]
        indices = P.compress(check, range(m))
        print "out-of-range xvec's = ", xvec[indices]
        print "out-of-range xvec indices = ", indices
        return False

    return True

#===============================================================================
def pchip_eval(x, y, m, xvec):

    """
    Evaluate the piecewise cubic Hermite interpolant with monoticity preserved
    
        x = array containing the x-data
        y = array containing the y-data
        m = slopes at each (x,y) point [computed to preserve monotonicity]
        xnew = new "x" value where the interpolation is desired
    
        x must be sorted low to high... (no repeats)
        y can have repeated values
    
    This works with either a scalar or vector of "xvec"
    """
    n = len(x)
    mm = len(xvec)

    if not x_is_okay(x, xvec):                                                  # Make sure there aren't problems with the input data
        print "pchip_eval2() - ill formed 'x' vector!!!!!!!!!!!!!"

        STOPpchip_eval2                                                         # Cause a hard crash...
    """
    Find the indices "k" such that x[k] < xvec < x[k+1]
    """
    
    xx = P.resize(x,(mm,n)).transpose()                                         # Create "copies" of "x" as rows in a mxn 2-dimensional vector
    xxx = xx > xvec

    z = xxx[:-1,:] - xxx[1:,:]                                                  # Compute column by column differences

    k = z.argmax(axis=0)                                                        # Collapse over rows...

    h = x[k+1] - x[k]                                                           # Create the Hermite coefficients
    t = (xvec - x[k]) / h

    h00 = (2 * t**3) - (3 * t**2) + 1                                           # Hermite basis functions
    h10 =      t**3  - (2 * t**2) + t
    h01 = (-2* t**3) + (3 * t**2)
    h11 =      t**3  -      t**2

    ynew = h00*y[k] + h10*h*m[k] + h01*y[k+1] + h11*h*m[k+1]                    # Compute the interpolated value of "y"

    return ynew

#===============================================================================
def pchip_init(x,y):
    
    """
    Evaluate the piecewise cubic Hermite interpolant with monoticity preserved
    
        x = array containing the x-data
        y = array containing the y-data
    
        x must be sorted low to high... (no repeats)
        y can have repeated values
    
        x input conditioning is assumed but not checked
    """
    
    n = len(x)

    delta = (y[1:] - y[:-1]) / (x[1:] - x[:-1])                                 # Compute the slopes of the secant lines between successive points

    m = P.zeros(n, dtype='d')                                                   # Initialize the tangents at every points as the average of the secants

    m[0] = delta[0]                                                             # At the endpoints - use one-sided differences
    m[n-1] = delta[-1]

    m[1:-1] = (delta[:-1] + delta[1:]) / 2.0                                    # In the middle - use the average of the secants

    """
    Special case: intervals where y[k] == y[k+1]

    Setting these slopes to zero guarantees the spline connecting
    these points will be flat which preserves monotonicity
    """
    
    indices_to_fix = P.compress((delta == 0.0), range(n))
    
#    print "zero slope indices to fix = ", indices_to_fix

    for ii in indices_to_fix:
        m[ii]   = 0.0
        m[ii+1] = 0.0

    alpha = m[:-1]/delta
    beta  = m[1:]/delta
    dist  = alpha**2 + beta**2
    tau   = 3.0 / P.sqrt(dist)

    """
    To prevent overshoot or undershoot, restrict the position vector
    (alpha, beta) to a circle of radius 3.  If (alpha**2 + beta**2)>9,
    then set m[k] = tau[k]alpha[k]delta[k] and m[k+1] = tau[k]beta[b]delta[k]
    where tau = 3/sqrt(alpha**2 + beta**2).
    """
    
    over = (dist > 9.0)                                                         # Find the indices that need adjustment
    indices_to_fix = P.compress(over, range(n))

#    print "overshoot indices to fix... = ", indices_to_fix

    for ii in indices_to_fix:
        m[ii]   = tau[ii] * alpha[ii] * delta[ii]
        m[ii+1] = tau[ii] * beta[ii]  * delta[ii]

    return m


#===============================================================================
def main():
    # Step function test...

    P.figure(2)
    P.title("pchip() step function test")

    x = exp_time                                                                # Create a step function (will demonstrate monotonicity)
    y = exp_conc

    xvec = P.arange(0.0, 720000.0)                                             # Interpolate using monotonic piecewise Hermite cubic spline
 
    m    = pchip_init(x,y)                                                      # Create the pchip slopes slopes

    yvec = pchip_eval(x, y, m, xvec)                                            # Interpolate...

    P.plot(x,    y,     'ro')                                                   # Plot the results
    P.plot(xvec, yvec,  'b')
    P.xlabel("Time (hours)")
    P.ylabel("[VOC] (ug/m3)")
    P.title("PCHIP")
    legends = ["Data", "pypchip()"]
    P.legend(legends, loc="upper right")

    P.show()
    #print x
    #print y    
    #print xvec
    #print yvec
    #print yvec3
    """
    Determination of Both Bounds for Initial Concentration. Must do one bound at
    a time:
        Co,u - add calculated tupper to end of exp_time
               add 0.0 to end of exp_conc
               change xvec to (0.0,tupper)
               run program and record Co,u
               
        Co,l - add calculated tlower to end of exp_time
               add 0.0 to end of exp_conc
               change xvec to (0.0,tlower)
               run program and record Co,l
    ================================================================================
    """
    Co = (air_flow_rate/(m_surf_area*m_thickness))*simps(yvec,xvec)
    print "Co = ",Co
main()
