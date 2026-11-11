import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# input variables

# max iterations for the solvers
maxIt = 10000
tol = 1e-8

# parameters for earth and sun
m_sun = 1.989e30
m_earth = 5.972e24
R = 1.496e8

# variable substitution
mu = m_earth / (m_earth + m_sun)

# initial guess
x = 0.5


# F(x) := Function for which the root is found
def F(x):
    return x ** 5 + (mu - 3) * x ** 4 + (3 - 2 * mu) * x ** 3 - mu * x ** 2 + 2 * mu * x - mu


# F(x)_prime := derivative of the function for which the root is found
def F_prime(x):
    return 5 * x ** 4 + 4 * (mu - 3) * x ** 3 + 3 * (3 - 2 * mu) * x ** 2 - 2 * mu * x + 2 * mu


x_iterations = []

# newton rhapson
cnt = 0
x_old = 0

while cnt < maxIt and np.abs(x - x_old) > tol:
    x_old = x

    # newton rhapson step
    x = x - F(x)/F_prime(x)

    # error calculation for later
    x_iterations.append(x)
    cnt += 1

newton_result = x

# calculate errors relative to
newton_errors = np.array(x_iterations) - newton_result

# quadr_errors := used to calculate the quadratic behavior of the error
newton_quadr_errors = np.zeros_like(newton_errors)
# initialise first element
newton_quadr_errors[0] = newton_errors[0]

# iterate n times to calculate the quadratic errors
for i in range(1, cnt):
    newton_quadr_errors[i] = newton_quadr_errors[i - 1] * newton_quadr_errors[i - 1]

# plot errors
plt.title("Error vs Iterations")
plt.xlabel("N")
plt.ylabel("Error")
plt.plot(np.linspace(0, len(newton_errors), len(newton_errors)), np.array(newton_errors), label="normal error")
plt.plot(np.linspace(0, len(newton_quadr_errors), len(newton_quadr_errors)), (np.array(newton_quadr_errors)), label="quadratic error")
plt.legend()
plt.show()

# bisection
# x element [0, 1]

# intervals := controls how many intervals the domain is split up into at each iteration
intervals = 10
start = 0
stop = 1

bracket_errors = []

cnt = 0

# loop while the error is larger than the treshhold
while stop - start > tol and cnt < maxIt:
    # split domain into subintervals
    intv = np.linspace(start, stop, intervals)

    # sign of first point
    minus = F(start) < 0

    stop_idx = 0

    # find the index for which the sign changes
    for i in range(1, len(intv)):
        # check sign change
        if (minus and (F(intv[i]) > 0)) or (not minus and (F(intv[i]) < 0)):
            stop_idx = i
            break

    # assign new upper and lower bound
    start = intv[stop_idx-1]
    stop = intv[stop_idx]

    # save error
    bracket_errors.append(stop - start)

    cnt += 1

bracketing_result = start

# create numpy array
bracket_errors = np.array(bracket_errors)

# plot the errors against each other
plt.title("Newton Rhapson vs Bracketing method convergence")
plt.xlabel("N")
plt.ylabel("Error")
plt.plot(np.linspace(0, len(bracket_errors), len(bracket_errors)), bracket_errors, label="bracketing")
plt.plot(np.linspace(0, len(newton_errors), len(newton_errors)), (np.array(newton_errors)), label="newton")
plt.legend()
plt.show()

print("##############################")
print(f"Root found for newton-rhapson: {newton_result}")
print(f"Root found for bracketing: {bracketing_result}")
print(f"Location of L1: {newton_result * R} km")
print("##############################")
