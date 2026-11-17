import matplotlib

matplotlib.use("TkAgg")

import numpy as np
import matplotlib.pyplot as plt

from gauss_legendre import gauss_legendre
from simpsons_rule import simpsons_rule
from trapezoidal_rule import trapezoidal_rule

N = 100
interval = np.linspace(np.pi, 0, N, endpoint=False)

T_cosh = np.zeros(N)
T_e = np.zeros(N)
T_cos = np.zeros(N)

# m = 1
# iterate over all N x datapoints
for i in range(N):
    a = interval[i]

    # create functions to integrate use the sub x = a sin(y)
    def f_cosh_trans(x):
        return a * np.cos(x) / np.sqrt(np.cosh(a) - np.cosh(a * np.sin(x)))


    def f_e_trans(x):
        return a * np.cos(x) / np.sqrt(np.exp(np.abs(a)) - np.exp(np.abs(a * np.sin(x))))


    def f_cos_trans(x):
        return a * np.cos(x) / np.sqrt(-np.cos(a) + np.cos(a * np.sin(x)))


    # integration
    T_cosh[i] = np.sqrt(8) * gauss_legendre(a=0, b=np.pi/2, n=4, f=f_cosh_trans)
    T_e[i] = np.sqrt(8) * gauss_legendre(a=0, b=np.pi/2, n=4, f=f_e_trans)
    T_cos[i] = np.sqrt(8) * gauss_legendre(a=0, b=np.pi/2, n=4, f=f_cos_trans)

# plotting
# plot functions
# plt.plot(interval, T_cosh)
# plt.title("V(x) = cosh(x)")
# plt.grid()
# plt.xlabel("a")
# plt.ylabel("T(a)")
# plt.show()

# plt.plot(interval, T_e)
# plt.title(r"V(x) = $e^{|x|}$")
# plt.xlabel("a")
# plt.grid()
# plt.ylabel("T(a)")
# plt.show()

plt.plot(interval, T_cos)
plt.title("V(x) = -cos(x)")
plt.xlabel("a")
plt.grid()
plt.ylabel("T(a)")
plt.show()
