import numpy as np


# a, b := Start- and Endpoint
# N := Number of Gridpoints
# f := function to evaluate
# endpoint := whether the last point should be included or not (for singularities)
def simpsons_rule(a, b, N, f, endpoint=True):
    # make sure N is dividable by 3, so that the interval can be split evenly
    N = N + N % 3 - 1

    # calculate grid spacing
    h = (b - a) / N

    # create evenly spaced gridpoints
    # included endpoint as it is a closed integration method
    interval = np.linspace(a, b, N, endpoint=endpoint)

    # summing variable for loop
    result = 0

    # calculate subintervals in loop, range from 0 to N - 3, so that last subinterval can be evaluated
    # make 3 steps at once
    for i in range(0, N - 3, 3):
        # add up subinterval according to formula
        result += 3 * h / 8 * (f(interval[i]) + 3 * f(interval[i + 1]) + 3 * f(interval[i + 2]) + f(interval[i + 3]))

    return result
