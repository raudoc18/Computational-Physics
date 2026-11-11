from ctypes import c_double, c_int, POINTER, byref, CDLL
import os

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# import c libraries
lib_dir = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(lib_dir, "find_roots.so")
c_lib = CDLL(lib_path)

# set parameter datatypes and return datatypes
c_lib.findRoots.argtypes = [c_double, c_double, c_double,
                             c_double, c_double, c_double, c_double, c_double, 
                             c_double, c_double, POINTER(c_int)]
c_lib.findRoots.restype = POINTER(c_double)

# set the number of iterations and the step (minimum distance between roots)
cnt = 0
N = 100000
step = 10e-4

# intervall
a = -10
b = 10

num_roots = []

# iterate N times
while cnt < N:

    # generate random coefficients
    args = np.random.rand(7)

    # root_cnt as a c integer (to pass it by reference)
    root_cnt = c_int(0)
    # call c function, returns a pointer of the double array of roots
    double_pointer = c_lib.findRoots(a, b, step, args[0], args[1], args[2], args[3], args[4], args[5], args[6], byref(root_cnt))

    # collect values
    root_cnt = root_cnt.value
    num_roots.append(root_cnt)

    cnt += 1

roots_arr = np.array(num_roots)

# make histogram to plot
counts, bins = np.histogram(roots_arr, bins=range(0, 7))

# plot the number of roots and their average
plt.bar(bins[:-1], counts, width=0.8, align="edge")
plt.xticks(bins)
plt.xlabel("Roots")
plt.ylabel("Frequency")
plt.title(f"N = {N}, Average = {float(np.average(roots_arr))}")
plt.show()