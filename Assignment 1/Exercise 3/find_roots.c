#include <stdlib.h>

// power function
double power(double x, int e) {
    double res = x;
    if (e == 0) {
        return 1;
    }
    for (int i = 1; i < e; i++) {
        res *= x;
    }
    return res;
}

// evaluates the polynomial at point x with the corresponding coefficients
double eval_poly(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6){
    return power(x, 6) * a6 + power(x, 5) * a5 + power(x, 4) * a4 + power(x, 3) * a3 + power(x, 2) * a2 + x * a1 + a0;
}

// findRoots := finds all roots in the intervall a, b by taking steps until x = b
// a,b := intervall
// a_n := coefficients of polynomial
// root_cnt := number of roots found
double* findRoots(double a, double b, double step, double a0, double a1, double a2, double a3, double a4, double a5, double a6, int* root_cnt){
    // eval first point
    double x = a;
    double current_sign = eval_poly(a, a0, a1, a2, a3, a4, a5, a6);

    // take step
    x += step;

    // initialise the roots array
    double *roots = (double*)malloc(6 * sizeof(double));

    // take steps until x is at b
    while (x <= b){
        // evaluate new point
        double new_sign = eval_poly(x, a0, a1, a2, a3, a4, a5, a6);
        // check if signs differ
        if ((current_sign > 0 && new_sign < 0) || (current_sign < 0 && new_sign > 0)){
            current_sign = new_sign;

            // count incremental
            roots[*root_cnt] = x;
            (*root_cnt)++;

            // if 6 roots have been found, stop as the polynomial can have 6 roots at max
            if (*root_cnt == 6){
                return roots;
            }
        }

        // take next step
        x += step;
    }

    return roots;
}