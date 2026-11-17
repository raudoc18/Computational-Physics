// power function
long double power(double x, int e) {
    double res = x;
    if (e == 0) {
        return 1;
    }
    for (int i = 1; i < e; i++) {
        res *= x;
    }
    return res;
}

// calculates nth legendre polynomial with recursive formula
double n_legendre(int n, double x) {
    if (n == 0){
        return 1;
    }
    else if (n == 1){
        return x;
    }

   return 1/((double)n) * ((2.0 * (double)n - 1.0) * x * n_legendre(n - 1, x) - ((double)n - 1.0) * n_legendre(n - 2, x));
}

// calculates the derivate of the nth legendre polynomial with recursive formula
double n_legendre_prime(int n, double x) {
    if (n == 0){
        return 0;
    }
    else if (n == 1){
        return 1;
    }

    double n_legendre_nm2 = n_legendre(n - 2, x);
    double n_legendre_nm1 = n_legendre(n - 1, x);
    double n_legendre_n = 1/((double)n) * ((2.0 * (double)n - 1.0) * x * n_legendre_nm1 - ((double)n - 1.0) * n_legendre_nm2);


    return 1/(power(x, 2) - 1) * ((double)n * x * n_legendre_n - (double)n * n_legendre_nm1);
}