#define functions_blinded_gm2_C
#include "functions_blinded_gm2.hpp"

double lhs_function_blinded_gm2_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}


double rhs_one_line(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    double r = P[0] * a2;
    return r;
}

double rhs_one_line_p1(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    double r = P[0] * a2 + 1.0;
    return r;
}

double rhs_line(int n, int Nvar, double* x, int Npar, double* P) {
    double a2 = x[0];
    double r = P[0] + P[1] * a2;
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
double rhs_a4TM(int n, int Nvar, double* x, int Npar, double* P) {

    double a2 = x[0];

    double r = P[0];
    switch (n) {
    case 0:
        r += P[1] * a2 + P[3] * a2 * a2;
        break;
    case 1:
        r += P[2] * a2;
        break;
    default:
        break;
    }

    return r;
}
double rhs_a4OS(int n, int Nvar, double* x, int Npar, double* P) {

    double a2 = x[0];

    double r = P[0];
    switch (n) {
    case 0:
        r += P[1] * a2 ;
        break;
    case 1:
        r += P[2] * a2 + P[3] * a2 * a2;
        break;
    default:
        break;
    }

    return r;
}
double rhs_a4(int n, int Nvar, double* x, int Npar, double* P) {

    double a2 = x[0];

    double r = P[0];
    switch (n) {
    case 0:
        r += P[1] * a2 + P[3] * a2 * a2;
        break;
    case 1:
        r += P[2] * a2 + P[4] * a2 * a2;
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
            X[j] += ZA[Njack - 1] * ZA[Njack - 1] * in[j][id_TM][t][0] - ZV[Njack - 1] * ZV[Njack - 1] * in[j][id_OS][t][0];
        }
    }

    return X;

}
double lhs_to_C80_to_Mphys(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[n]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j];

    return r;
}

double lhs_TM_m_OS(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[0]][j] - gjack.en[e].jack[fit_info.corr_id[1]][j];

    return r;
}

double lhs_TM_over_OSto_C80_to_Mphys(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[0]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j];
    r /= gjack.en[e].jack[fit_info.corr_id[1]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j];

    return r;
}

double lhs_aveTMOS_to_C80_to_Mphys(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[0]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j];
    r += gjack.en[e].jack[fit_info.corr_id[1]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j];
    r /= 2;

    return r;
}