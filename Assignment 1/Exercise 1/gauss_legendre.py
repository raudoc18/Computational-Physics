import ctypes
import os
import numpy as np

# include c code into python for faster calculation of legendre polynomials
lib_dir = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(lib_dir, "lib_legendre_poly.so")
c_lib = ctypes.CDLL(lib_path)

c_lib.n_legendre.argtypes = [ctypes.c_int, ctypes.c_double]
c_lib.n_legendre.restype = ctypes.c_double

c_lib.n_legendre_prime.argtypes = [ctypes.c_int, ctypes.c_double]
c_lib.n_legendre_prime.restype = ctypes.c_double

c_lib.fact.argtypes = [ctypes.c_int]
c_lib.fact.restype = ctypes.c_int
c_lib.binom.argtypes = [ctypes.c_int, ctypes.c_int]
c_lib.binom.restype = ctypes.c_int
c_lib.power.argtypes = [ctypes.c_double, ctypes.c_int]
c_lib.power.restype = ctypes.c_double

# transform of function to legendre space
# a, b := start- and endpoint
# f := function to be transformed
def F(a, b, x, f):
    return (b - a) / 2 * f((b - a) / 2 * x + (b + a) / 2)


# a, b := start- and endpoint
# n := degree of polynomials
# f := integrand function
# maxIt := maximum number of iterations to find root
# tol := error tolerance for root
def gauss_legendre(a, b, n, f, maxIt=10000, tol=1e-6):
    integral_res = 0

    for i in range(1, n + 1):
        # find ith root of legendre polynomial of degree n
        # approximate ith root
        x = np.cos(np.pi * (4 * i - 1) / (4 * n + 2))

        # initial old value to 2 because we are between -1 and 1
        xold = 2

        # iteration count
        cnt = 0

        # newton rhapson to find ith root
        while np.abs(xold - x) > tol and cnt < maxIt:
            xold = x
            x = x - c_lib.n_legendre(n, x) / c_lib.n_legendre_prime(n, x)
            cnt += 1

        # calculate weight i according to formula
        w_i = 2 / ((1 - x ** 2) * (c_lib.n_legendre_prime(n, x) ** 2))

        # calculate F at position xi
        F_i = F(a, b, x, f)

        integral_res += w_i * F_i

    return integral_res

