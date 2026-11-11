import numpy as np


# a, b := Start- and Endpoint
# N := Number of Gridpoints
# f := function to evaluate
# endpoint := whether the last point should be included or not (for singularities)
def trapezoidal_rule(a, b, N, f, endpoint=True):
    # create evenly spaced gridpoints
    # included endpoint as it is a closed integration method
    interval = np.linspace(a, b, N, endpoint=endpoint)

    # calculate the grid spacing
    h = (b - a) / N

    # calculate the sum of trapezoids, use array slicing for right part
    result = h / 2 * (f(interval[0]) + f(interval[-1])) + h * np.sum(f(interval[1:-1]))

    return result
