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

double* unblind_combo(double**** in, double* ZA, double* ZV, struct fit_type fit_info) {
    int Njack = fit_info.Njack;
    int Tmax = fit_info.T > 31 ? 31 : fit_info.T;
    int id_TM = fit_info.corr_id[0];
    int id_OS = fit_info.corr_id[1];

    double* X = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        X[j] = 0;
        for (int t = 1;t < Tmax;t++) {
            X[j] += ZA[Njack-1] * ZA[Njack-1] * in[j][id_TM][t][0] - ZV[Njack-1] * ZV[Njack-1] * in[j][id_OS][t][0];
        }
    }

    return X;

}
