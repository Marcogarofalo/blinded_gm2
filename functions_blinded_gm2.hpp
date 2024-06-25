#ifndef functions_blinded_gm2_H
#define functions_blinded_gm2_H
#include "non_linear_fit.hpp"

double lhs_function_blinded_gm2_eg(int j, double**** in, int t, struct fit_type fit_info);
double rhs_two_line(int n, int Nvar, double* x, int Npar, double* P);
double* unblind_combo(double**** in, double* ZA, double* ZV, struct fit_type fit_info);
#endif
