#define functions_blinded_gm2_C
#include "functions_blinded_gm2.hpp"

double lhs_function_blinded_gm2_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}



double rhs_two_line(int n, int Nvar, double* x, int Npar, double* P) {

    double a2 = x[0];

    double r = P[0];
    switch (n) {
    case 0:
        r += P[1] * a2;
        break;
    case 1:
        r += P[2] * a2;
        break;
    default:
        break;
    }

    return r;
}

