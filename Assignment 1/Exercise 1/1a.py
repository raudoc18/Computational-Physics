import matplotlib
matplotlib.use("TkAgg")

import numpy as np
import matplotlib.pyplot as plt

from gauss_legendre import gauss_legendre
from simpsons_rule import simpsons_rule
from trapezoidal_rule import trapezoidal_rule

# constants
a = 0
b = 1
# define number of Ns and ns
N = np.logspace(7, 17, base=2, num=10, dtype=int)
n = np.logspace(1, 5, num=5, base=2, dtype=int)
print(n)


# define integrand
def f(x):
    return 1 / np.sqrt(1 - x ** 4)


# define transformed integrand
a_trans = 0
b_trans = np.pi / 2


def f_trans(x):
    return 1 / np.sqrt(1 + np.sin(x) ** 2)


# create arrays to compare results afterward
results_trapezoidal = []
results_simpsons = []
results_trapezoidal_trans = []
results_simpsons_trans = []
results_gauss = []
results_gauss_trans = []

for Ni in N:
    # append to result arrays in each iteration
    # exclude endpoint to avoid divergence
    results_trapezoidal.append(trapezoidal_rule(a, b, Ni, f, False))
    results_simpsons.append(simpsons_rule(a, b, Ni, f, False))

    results_trapezoidal_trans.append(trapezoidal_rule(a_trans, b_trans, Ni, f_trans, True))
    results_simpsons_trans.append(simpsons_rule(a_trans, b_trans, Ni, f_trans, True))

for ni in n:
    # append to result arrays in each iteration
    results_gauss.append(gauss_legendre(a, b, ni, f))
    results_gauss_trans.append(gauss_legendre(a_trans, b_trans, ni, f_trans))

# create numpy arrays to make calculations possible
results_trapezoidal = np.array(results_trapezoidal)
results_trapezoidal_trans = np.array(results_trapezoidal_trans)

results_simpsons = np.array(results_simpsons)
results_simpsons_trans = np.array(results_simpsons_trans)

results_gauss = np.array(results_gauss)
results_gauss_trans = np.array(results_gauss_trans)

# plot errors according to error formula
# plt.plot(N[1:], np.abs(results_trapezoidal[1:] - results_trapezoidal[:-1]) / 3, label="normal integral")
# plt.plot(N[1:], np.abs(results_trapezoidal_trans[1:] - results_trapezoidal_trans[:-1]) / 3,
#          label="transformed integral")
# plt.legend()
# plt.xlabel("N")
# plt.ylabel("Error")
# plt.title("Trapezoidal Error")
#
# plt.show()
#
# # plot errors according to error formula
# plt.plot(N[1:], np.abs(results_simpsons[1:] - results_simpsons[:-1]) / 15, label="normal integral")
# plt.plot(N[1:], np.abs(results_simpsons_trans[1:] - results_simpsons_trans[:-1]) / 15, label="transformed integral")
# plt.legend()
# plt.xlabel("N")
# plt.ylabel("Error")
# plt.title("Simpsons Error")
#
# plt.show()
#
# # plot errors according to error formula
# plt.plot(n[1:], np.abs(results_gauss[1:] - results_gauss[:-1]), label="normal integral")
# plt.plot(n[1:], np.abs(results_gauss_trans[1:] - results_gauss_trans[:-1]), label="transformed integral")
# plt.legend()
# plt.xlabel("n")
# plt.ylabel("Error")
# plt.title("Gauss-Legendre Error")
#
# plt.show()

print("#############################")
print("Results of integration: ")
print(f"Trapezoidal: T(1) = {results_trapezoidal_trans[-1]}, N = {N[-1]}")
print(f"Simpsons: T(1) = {results_simpsons_trans[-1]}, N = {N[-1]}")
print(f"Trapezoidal: T(1) = {results_gauss_trans[-1]}, n = {n[-1]}")
print("#############################")
