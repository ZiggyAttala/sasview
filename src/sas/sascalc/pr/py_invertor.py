

import numpy as np
import sys
import math
import time
import copy
import os
import re
import logging
"""
Implementing C invertor functions in Python
"""
def ortho_transformed_py(d_max, n, q):
    return 8.0 * math.pow(math.pi, 2.0)/q * d_max * n * math.pow(-1.0, n + 1) * math.sin(q*d_max) / ( math.pow(math.pi*n, 2.0)
    - math.pow(q * d_max, 2.0) )

def ortho_transformed_smeared_py(d_max, n, height, width, q, npts):
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
    
    print("height: ", height)
    print("width: ", width)
    print("npts", npts)
    n_height = (1, npts)[height>0]
    n_width = (1, npts)[width>0]
    
    print("n_height", n_height)
    print("n_width", n_width)

    count_w = 0.0
    print(fnpts)

    z_func = lambda j,height=height,fnpts=fnpts : (0.0, height/fnpts*float(j))[height>0]
    
    j_iter = np.arange(n_height)
    i_iter = np.arange(n_width)

    z_result = map(z_func, j_iter)
    
    z_values = np.fromiter((z_result), dtype = np.float)
    print("z values: ", z_values)
    y_func = lambda i, width=width, fnpts=fnpts : (0.0, -width/2.0+width/fnpts*float(i))[width>0]
    
    y_result = map(y_func, i_iter)
    y_values = np.fromiter((y_result), dtype = np.float)

    print(z_values)
    print(y_values)

    calc_func = lambda y,z,q=q : (((q - y) * (q - y) + z * z))

    calc_result = map(calc_func, y_values, z_values)

    calc_values = np.fromiter((calc_result), dtype = np.float)
    print(calc_values)
    result_to_transform = filter(lambda x : x > 0, calc_values)

    vals_to_transform = np.fromiter((result_to_transform), dtype = np.float)

    
    transformed_values = map(lambda x : ortho_transformed_py(d_max, n, math.sqrt(x)), vals_to_transform)
    
    count_w = calc_values.shape[0]
    sum = np.sum(np.fromiter(vals_to_transform, dtype = np.float))
    print(sum)
    print(count_w)
    return sum/count_w

def ortho_transformed_smeared(d_max, n, height, width, q, npts):
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

    n_height = (1, npts)[height>0]
    n_width = (1, npts)[width>0]

    count_w = 0.0

    for j in range(0,n_height):
        if(height>0):
            z = height/fnpts* float(j)
        else:
            z = 0.0
        for i in range(0, n_width):
            if(width>0):
                y = -width/2.0+width/fnpts* float(i)
            else:
                y = 0.0
            if (((q - y) * (q - y) + z * z) > 0.0):
                count_w += 1.0
                sum += ortho_transformed_py(d_max, n, math.sqrt((q - y)*(q - y) + z * z))
    print(sum)
    print(count_w)
    return sum / count_w


start_time = time.time()
print("Result normal: ", ortho_transformed_smeared(1, 1, 1000, 1000, 1, 5))
end_time = time.time()
print("Time 1:", end_time - start_time)

start_time_py = time.time()
print("Result py:", ortho_transformed_smeared_py(1, 1, 1000, 1000, 1, 5))
end_time_py = time.time()
print("Time 2: ", end_time_py - start_time_py)

