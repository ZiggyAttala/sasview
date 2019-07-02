"""
Starting to convert invertor.c, eventually the Cinvertor class to python

At the moment experimenting with different methods of implementing individual
functions.
"""
import numpy as np
import sys
import math
import time
import copy
import os
import re
import logging
import time
from numba import jit, njit, vectorize, float64, guvectorize

#class stub for final Pinvertor class
#taking signatures from Cinvertor.c and docstrings
class Pinvertor:
    def __init__(self):
        pass
    def residuals(self, args):
        """
        Function to call to evaluate the residuals\n
	    for P(r) inversion\n
	    @param args: input parameters\n
	    @return: list of residuals
        """
        pass

#Private Methods


@njit
def pr_sphere(R, r):
    """
    P(r) of a sphere, for test purposes

    @param R: radius of the sphere
    @param r: distance, in same units as the radius
    @return: P(r)
    """
    if(r <= 2.0*R):
        return 12.0* (0.5*r/R)**2 * (1.0-0.5*r/R)**2 * (2.0+0.5*r/R)
    else:
        return 0

@njit
def ortho(d_max, n, r):
    """
    Orthogonal Functions:
    B(r) = 2r sin(pi*nr/d)
    """
    return 2.0 * r * np.sin(np.pi * n * r/d_max)

@njit()
def ortho_transformed(d_max, n, q):
    """
    Fourier transform of the nth orthagonal function
    
    With vectorize time was 
    \@vectorize() ~= 1.4e-05
    \@njit() ~= 3e-06
    \@njit(parallel=True) - Compiler returns no transformation for parallel execution possible
    and time was the same as @njit()
    """
    #return 8.0*(np.pi)**2/q * d_max * n * (-1.0)**(n+1) * math.sin(q*d_max) / ( (math.pi*n)**2 - (q*d_max)**2 )
    pi = np.pi
    qd = q * (d_max/pi)
    return ( 8.0 * d_max**2/pi * n * (-1.0)**(n+1) ) * np.sinc(qd) / (n**2 - qd**2)



@njit()
def ortho_transformed_smeared(d_max, n, height, width, q, npts):
    """
    Slit-smeared Fourier transform of the nth orthagonal function.
    Smearing follows Lake, Acta Cryst. (1967) 23, 191.

    \@njit() - 5 - time roughly 4.5e-05
    \@njit() - npts 1000 - 0.031
    \@njit(parallel = True) - no transformation possible same time.
    \@vectorize([float64(float64, float64, float64, float64, float64, float64)])  npts = 5 ~= 1.6e-05
    npts = 1000 ~= 0.02
    """
    y = 0
    z = 0
    sum = 0

    i = 0
    j = 0
    n_height = 0
    n_width = 0
    count_w = 0
    fnpts = 0
    sum = 0.0
    fnpts = float(npts-1.0)
    if(height > 0):
        n_height = npts
    else:
        n_height = 1
    if(width > 0):
        n_width = npts
    else:
        n_width = 1

    count_w = 0.0

    #Pre compute dz, y0 and dy
    dz = height/(npts-1)
    y0 = -0.5*width
    dy = width/(npts-1)

    for j in range(0, n_height):
        zsq = (j * dz)**2

        for i in range(0, n_width):
            y = y0 + i*dy
            qsq = (q - y)**2 + zsq
            count_w += qsq > 0
            sum += ortho_transformed(d_max, n, np.sqrt(qsq))
    return sum / count_w

@njit()
def ortho_derived(d_max, n, r):
    """
    First derivative in of the orthogonal function dB(r)/dr
    """
    pinr = pi * n * r/d_max
    return 2.0 * np.sin(pinr) + 2.0 * r * np.cos(pinr)

@njit()
def iq(pars, d_max, n_c, q):
    """
    Scattering intensity calculated from the expansion

    basic python ~=0.00014 with array of size 20
    vectorised operations ~= 0.0002
    using njit, no parallel possible ~= 1.3e-05
    """
    sum = 0.0
    for i in range(0, n_c):
        sum += pars[i] * ortho_transformed(d_max, i + 1, q)

    return sum

@jit()
def iq_smeared_vec(pars, d_max, n_c, height, width, q, npts):
    """
    Scattering intensity calculated from expansion,
    slit smeared.
    for test data of size 20 ~= 0.12
    """
    sum = 0.0

    @vectorize("float64(float64, float64)", nopython = True)
    def f(x, y):
        return x * ortho_transformed_smeared(d_max, y + 1, height, width, q, npts)

    i = np.arange(n_c)

    result = f(pars, i)
    return np.sum(result)
    #for i in range(0, n_c):
     #   sum += pars[i] * ortho_transformed_smeared(d_max, i + 1, height, width, q, npts)
    #return sum

@njit()
def iq_smeared_njit(pars, d_max, n_c, height, width, q, npts):
    """
    Scattering intensity calculated from expansion,
    slit smeared.
    test data 20, 0.0003
    """
    sum = 0.0

    for i in range(0, n_c):
        sum += pars[i] * ortho_transformed_smeared(d_max, i + 1, height, width, q, npts)
    return sum

def iq_smeared_basic(pars, d_max, n_c, height, width, q, npts):
    """
    Scattering intensity calculated from expansion,
    slit smeared.
    test data 20, 0.0005
    """
    sum = 0.0

    for i in range(0, n_c):
        sum += pars[i] * ortho_transformed_smeared(d_max, i + 1, height, width, q, npts)
    return sum

@njit()
def pr(pars, d_max, n_c, r):
    """
    P(r) calculated from the expansion
    """
    sum = 0.0
    for i in range(0, n_c):
        sum += pars[i] * ortho(d_max, i+1, r)
    return sum




#testing

def demo_ot():
    print(ortho_transformed(1,1,0))

def demo():
    tests = 10
    for i in range(0, tests): 
        start = time.clock()
        x = iq_smeared_vec(np.arange(20), 3, 20, 10, 10, 20, 20)
        end = time.clock()
        print(x)
        print("Time elapsed py : %s" % (end - start))

if(__name__ == "__main__"):
    demo()