import numpy as np

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from gauss_legendre import gauss_legendre

# eta, used to scale the steps taken
eta = 0.0001

# integrand_a function => returns the function of the integrand with specified alpha values
def integrand_a(alpha):
    def I(x):
        return np.cos(x) / np.sqrt(1 - np.sin(x) ** alpha)

    return I


# integrand_a_prime function => returns the function of the integrand of the derivative with respect to alpha
# with specified alpha values
def integrand_a_prime(alpha):
    def I(x):
        return 0.5 * np.cos(x) * (1 - np.sin(x)**alpha)**(-3/2) * (np.log(2)/2**alpha + np.sin(x)**alpha * np.log(np.sin(x)))

    return I

# integrand_a_half function => returns the function of the integrand of the a/2 integral with specified alpha values
def integrand_a_half(alpha):
    def I(x):
        return np.cos(x) / np.sqrt(0.5**alpha - np.sin(x)**alpha)

    return I

# integrand_a_half_prime function => returns the function of the integrand of the a/2 integral of the derivative with respect to alpha
# with specified alpha values
def integrand_a_half_prime(alpha):
    def I(x):
        return 0.5 * np.cos(x) * ((0.5**alpha - np.sin(x)**alpha)**(-3/2) * (np.sin(x)**alpha * np.log(np.sin(x))))

    return I


def gradient_decent():
    # choose a starting value for alpha
    alpha = 1

    # tolerance for the alpha error
    tol = 1e-3

    # initial value for grad_alpha to get into the loop
    grad_alpha = 1

    # repeat until grad_alpha is really low
    while np.linalg.norm([grad_alpha]) > tol:

        # calculate integrals using gauss legendre
        I_a = gauss_legendre(0, np.pi/2, n=4, f=integrand_a(alpha=alpha))
        I_a_prime = gauss_legendre(0, np.pi / 2, n=4, f=integrand_a_prime(alpha=alpha))
        I_a_half = gauss_legendre(0, np.arcsin(0.5), n=4, f=integrand_a_half(alpha=alpha))
        I_a_half_prime = gauss_legendre(0, np.arcsin(0.5), n=4, f=integrand_a_half_prime(alpha=alpha))

        # calculate the derivative with respect to alpha
        grad_alpha = (I_a_prime * I_a_half - I_a * I_a_half_prime)/I_a_half ** 2

        # make a step in grad_alpha direction
        alpha = alpha + grad_alpha * eta
        print(alpha)

    return alpha

alpha = gradient_decent()

# plot results

x = np.linspace(0, 1, 1000)
plt.plot(x, np.abs(x**alpha))
plt.title(r"Minimizing $\frac{T(a)}{T(a/2)} through V(x)$")
plt.grid()
plt.xlabel("x")
plt.ylabel("V(x)")
plt.show()
