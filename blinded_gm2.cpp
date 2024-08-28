#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"
#include "functions_blinded_gm2.hpp"

#include "functions_amu.hpp"
constexpr double Mrho_MeV = 775.0;
// constexpr double Mrho_MeV = 770.0;
constexpr double Mrho_MeV_err = 0.4;

// constexpr double grhopipi = 6.06;
// constexpr double grhopipi_err = 0.03;

constexpr double grhopipi = 5.95;
// constexpr double grhopipi = 5.5;
constexpr double grhopipi_err = 0.0001;


struct kinematic kinematic_2pt;

void weighted_average(double* out, double* in1, double* in2) {
    double w1 = 1.0 / myres->comp_error(in1);
    double w2 = 1.0 / myres->comp_error(in2);
    double se = w1 + w2;
    w1 /= se;
    w2 /= se;
    myres->linear_comb(out, w1, in1, w2, in2);
}

generic_header read_head(FILE* stream) {
    generic_header header;
    size_t Nconfs, Nhits, TT, Nsubs;
    int i = 0;
    i += fread(&Nconfs, sizeof(size_t), 1, stream);
    i += fread(&Nhits, sizeof(size_t), 1, stream);
    i += fread(&TT, sizeof(size_t), 1, stream);
    i += fread(&Nsubs, sizeof(size_t), 1, stream);
    std::cout << "Nconfs: " << Nconfs << std::endl;
    std::cout << "TT: " << TT << " " << TT / 2 + 1 << std::endl;
    std::cout << "Nhits: " << Nhits << std::endl;
    std::cout << "Nsubs: " << Nsubs << std::endl;
    error(i == 6, 1, "read_head", " counter %d, different from 6");

    int here = ftell(stream);
    fseek(stream, 0, SEEK_END);
    int size = ftell(stream);
    printf("Nconfs found: %ld\n", (size - here) / ((TT / 2 + 1) * sizeof(double)));
    fseek(stream, here, SEEK_SET);
    header.Nconf = Nconfs * Nsubs;
    error(header.Nconf != (size - here) / ((TT / 2 + 1) * sizeof(double)), 1, "read_head", "Nconf not match size file ");
    header.Njack = Nconfs;
    header.T = TT;
    header.L = TT / 2;
    return header;
}

void check_head(generic_header header, FILE* stream) {

    size_t Nconfs, Nhits, TT, Nsubs;
    int i = 0;
    i += fread(&Nconfs, sizeof(size_t), 1, stream);
    i += fread(&Nhits, sizeof(size_t), 1, stream);
    i += fread(&TT, sizeof(size_t), 1, stream);
    i += fread(&Nsubs, sizeof(size_t), 1, stream);
    std::cout << "Nconfs: " << Nconfs << std::endl;
    std::cout << "TT: " << TT << " " << TT / 2 + 1 << std::endl;
    std::cout << "Nhits: " << Nhits << std::endl;
    std::cout << "Nsubs: " << Nsubs << std::endl;
    error(header.Nconf != Nconfs * Nsubs, 1, "check_head", " Nconf * Nsub does not match read %d stored %ld", Nconfs, header.Nconf);
    error(header.T != TT, 1, "check_head", " T does not match read %d stored %ld", TT, header.T);
    error(header.Njack != Nconfs, 1, "check_head", " Nconf does not match read %d stored %d", Nconfs / Nsubs, header.Njack);
    error(i == 6, 1, "read_head", " counter %d, different from 6");

}
void write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[2]);         // path
    mysprintf(option[4], NAMESIZE, argv[6]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, argv[3]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[0];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[0];
    file_head.k[3] = head.mus[0];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

void read_twopt(FILE* stream, double*** to_write, generic_header head, int TMOS, int Nsub) {
    // write your function to read the data
    int fi = 0;
    // printf("aaaaaaaaaaaaaaaaa %d %d %d\n", head.Nconf / head.Njack, head.Nconf, head.Njack);
    for (int t = 0; t < head.T / 2 + 1;t++) {
        for (int k = 0; k < Nsub; k++) {
            double a, b;
            // fi += fscanf(stream, "%lf  %lf", &a, &b);
            fi += fread(&a, sizeof(double), 1, stream);
            to_write[TMOS][t][0] += a;
            // to_write[TMOS][t][1] += 0;
            // printf("%g\n", to_write[TMOS][t][0]);
        }
        to_write[TMOS][t][0] /= (double)Nsub;
        // to_write[TMOS][t][1] /= (double)Nsub;
    }
    error(fi != Nsub * (head.T / 2 + 1), 1, "bin_intoN", "invalid counter %d ", fi);
    for (int t = 1; t < head.T / 2 + 1;t++) {
        to_write[TMOS][head.T - t][0] /= to_write[TMOS][t][0];
        // to_write[TMOS][head.T - t][1] /= to_write[TMOS][t][1];
    }
}

int main(int argc, char** argv) {
    error(argc != 8, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac -p path file -bin $bin  jack/boot  mu\n "
        "separate path and file please");

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[6]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    mysprintf(namefile, NAMESIZE, "%s/%s/VkVk_tm.bin", option[3], option[6]);
    FILE* infile = open_file(namefile, "r");

    mysprintf(namefile, NAMESIZE, "%s/%s/VkVk_OS.bin", option[3], option[6]);
    FILE* infile_OS = open_file(namefile, "r");
    //////////////////////////////////// read and setup header
    generic_header head = read_head(infile);
    check_head(head, infile_OS);
    int Nconf, Nsub, Nskip;
    if (strcmp(option[6], "C80") == 0) {
        Nconf = head.Nconf - 218;
        Nsub = head.Nconf / head.Njack;
        Nskip = 218;
        head.Njack -= 218;
    }
    else {
        Nconf = head.Nconf;
        Nsub = head.Nconf / head.Njack;
    }
    // we need to init mu before pass the info to the global header
    head.mus.resize(1);
    head.mus[0] = std::atof(argv[7]);
    init_global_head(head);

    head.print_header();

    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[5]);
    int Neff = bin;
    error(bin < 10, 1, "main", "we are binning into Nbin, set a number >10");
    int Njack;
    if (strcmp(argv[6], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[6], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[7]= %s is not jack or boot", argv[7]);
    }
    // now Njack need to be the number of jacks
    head.Nconf = head.Njack;
    head.Njack = Njack;
    head.ncorr = 1;
    head.size = head.ncorr * head.T * 2;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4],
        option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    head.ncorr = 4;
    //////////////////////////////////// confs
    double**** data = calloc_corr(confs, head.ncorr, head.T);
    double**** data_1 = calloc_corr(218, head.ncorr, head.T);
    double**** data_bin_1;
    double**** conf_jack_1;

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    if (strcmp(option[6], "C80") == 0) {
        for (int iconf = 0; iconf < 218; iconf++) {
            read_twopt(infile, data_1[iconf], head, 0, Nsub);
            read_twopt(infile_OS, data_1[iconf], head, 1, Nsub);
        }
        data_bin_1 = bin_intoN(data_1, head.ncorr, head.T, 218, bin);
        // (confs, head.ncorr, head.T, data, bin);
        free_corr(218, head.ncorr, head.T, data_1);
        conf_jack_1 = myres->create(Neff, head.ncorr, head.T, data_bin_1);
        free_corr(Neff, head.ncorr, head.T, data_bin_1);
    }
    for (int iconf = 0; iconf < confs; iconf++) {
        read_twopt(infile, data[iconf], head, 0, Nsub);
        read_twopt(infile_OS, data[iconf], head, 1, Nsub);
    }


    double**** data_bin = bin_intoN(data, head.ncorr, head.T, confs, bin);//binning(confs, head.ncorr, head.T, data, bin);
    free_corr(confs, head.ncorr, head.T, data);
    /////////////////////////////////// read the three corr
    mysprintf(namefile, NAMESIZE, "/home/garofalo/analysis/g-2_new_stat/Vkvk_cont/%d_m%.5f/OPPOR", head.L, head.mus[0]);
    printf("reading: %s   \n", namefile);
    FILE* OPPOR = open_file(namefile, "r+");
    int ir = 0;
    ir += fscanf(OPPOR, "%*[^\n]\n");// skip a line
    for (int t = 0; t < head.T / 2 + 1;t++) {
        int tr;
        double a, b;
        ir += fscanf(OPPOR, "%d   %lf  %lf  %lf\n", &tr, &data_bin[0][2][t][0], &a, &b);
        error(tr != t, 1, "read free corr", "expected t=%d   , read t=%d", tr, t);
        for (int j = 1;j < bin;j++) {
            data_bin[j][2][t][0] = data_bin[0][2][t][0];
        }
    }
    mysprintf(namefile, NAMESIZE, "/home/garofalo/analysis/g-2_new_stat/Vkvk_cont/%d_m%.5f/SAMER", head.L, head.mus[0]);
    printf("reading: %s   \n", namefile);
    FILE* SAMER = open_file(namefile, "r+");
    ir += fscanf(SAMER, "%*[^\n]\n");// skip a line
    for (int t = 0; t < head.T / 2 + 1;t++) {
        int tr;
        double a, b;
        ir += fscanf(SAMER, "%d   %lf  %lf  %lf\n", &tr, &data_bin[0][3][t][0], &a, &b);
        error(tr != t, 1, "read free corr", "expected t=%d   , read t=%d", tr, t);
        for (int j = 1;j < bin;j++) {
            data_bin[j][3][t][0] = data_bin[0][3][t][0];
        }
    }
    double**** conf_jack = myres->create(Neff, head.ncorr, head.T, data_bin);
    free_corr(Neff, head.ncorr, head.T, data_bin);
    //////////////////////////////////////////////////////////////
    // read a,ZA,ZV
    //////////////////////////////////////////////////////////////
    double mean, err;
    int seed;
    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double* a_fm = myres->create_fake(mean, err, seed);
    double* a_MeV = myres->create_copy(a_fm);
    myres->div(a_MeV, a_fm, hbarc);
    line_read_param(option, "ZA", mean, err, seed, namefile_plateaux);
    double* ZA = myres->create_fake(mean, err, seed);
    line_read_param(option, "ZV", mean, err, seed, namefile_plateaux);
    double* ZV = myres->create_fake(mean, err, seed);
    line_read_param(option, "Mpi", mean, err, seed, namefile_plateaux);
    double* Mpi = myres->create_fake(mean, err, seed);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
        option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
        option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
        option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
        option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3],
        option[6]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3],
        option[6]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");


    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;

    for (int icorr = 0; icorr < head.ncorr; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;
    constexpr double q2ud = 5.0 / 9.0;
    double (*int_scheme)(int, int, double*);
    int_scheme = integrate_simpson38;
    int isub = -2;

    // 
    write_jack(a_fm, Njack, jack_file);
    write_jack(ZA, Njack, jack_file);
    write_jack(ZV, Njack, jack_file);
    write_jack(Mpi, Njack, jack_file);
    //////////////////////////////////////////////////////////////
    // no bounding
    //////////////////////////////////////////////////////////////
    double* amu_TM = compute_amu_full(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{full}_(TM)", resampling, isub);
    write_jack(amu_TM, Njack, jack_file);
    printf("amu_TM = %g  %g \n", amu_TM[Njack - 1], myres->comp_error(amu_TM));
    double* amu_OS = compute_amu_full(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{full}_(OS)", resampling, isub);
    write_jack(amu_OS, Njack, jack_file);
    printf("amu_OS = %g  %g \n", amu_OS[Njack - 1], myres->comp_error(amu_OS));

    //////////////////////////////////////////////////////////////
    // bounding BMW
    //////////////////////////////////////////////////////////////
    int ibound = 2;

    double* amu_TM_bound = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound}_(TM)",
        resampling, isub, ibound, Mpi, nullptr, head.L);
    write_jack(amu_TM_bound, Njack, jack_file);
    printf("amu_TM(bound) = %g  %g \n", amu_TM_bound[Njack - 1], myres->comp_error(amu_TM_bound));

    double* amu_OS_bound = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound}_(OS)",
        resampling, isub, ibound, Mpi, nullptr, head.L);
    write_jack(amu_OS_bound, Njack, jack_file);
    printf("amu_OS(bound) = %g  %g \n", amu_OS_bound[Njack - 1], myres->comp_error(amu_OS_bound));

    //////////////////////////////////////////////////////////////
    // bounding meff_t
    //////////////////////////////////////////////////////////////
    ibound = 3;

    double* amu_TM_bound_meff_t = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff_t}_(TM)",
        resampling, isub, ibound, Mpi, nullptr, head.L);
    write_jack(amu_TM_bound_meff_t, Njack, jack_file); check_correlatro_counter(8);
    printf("amu_TM(bound_meff_t) = %g  %g \n", amu_TM_bound_meff_t[Njack - 1], myres->comp_error(amu_TM_bound_meff_t));

    double* amu_OS_bound_meff_t = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff_t}_(OS)",
        resampling, isub, ibound, Mpi, nullptr, head.L);
    write_jack(amu_OS_bound_meff_t, Njack, jack_file); check_correlatro_counter(9);
    printf("amu_OS(bound_meff_t) = %g  %g \n", amu_OS_bound_meff_t[Njack - 1], myres->comp_error(amu_OS_bound_meff_t));

    //////////////////////////////////////////////////////////////
    // bounding meff
    //////////////////////////////////////////////////////////////
    ibound = 1;

    double* M_VKVK_TM = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 0, "M_VKVK_TM", M_eff_T, jack_file);

    double* M_VKVK_OS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 1, "M_VKVK_OS", M_eff_T, jack_file);

    double* amu_TM_bound_meff = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff}_(TM)",
        resampling, isub, ibound, Mpi, M_VKVK_TM, head.L);
    write_jack(amu_TM_bound_meff, Njack, jack_file);
    printf("amu_TM(bound_meff) = %g  %g \n", amu_TM_bound_meff[Njack - 1], myres->comp_error(amu_TM_bound_meff));

    double* amu_OS_bound_meff = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff}_(OS)",
        resampling, isub, ibound, Mpi, M_VKVK_OS, head.L);
    write_jack(amu_OS_bound_meff, Njack, jack_file);
    printf("amu_OS(bound_meff) = %g  %g \n", amu_OS_bound_meff[Njack - 1], myres->comp_error(amu_OS_bound_meff));

    double* amu_TM_bound_1;
    double* amu_OS_bound_1;
    double* amu_TM_bound_meff_t_1;
    double* amu_OS_bound_meff_t_1;
    double* amu_TM_bound_meff_1;
    double* amu_OS_bound_meff_1;
    double* zeros = (double*)calloc(Njack, sizeof(double));

    if (strcmp(option[6], "C80") == 0) {
        //////////////////////////////////////////////////////////////
        // bounding BMW
        //////////////////////////////////////////////////////////////
        int ibound = 2;

        amu_TM_bound_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound}_(TM)_1",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_TM_bound_1, Njack, jack_file);
        printf("amu_TM(bound)_1 = %g  %g \n", amu_TM_bound_1[Njack - 1], myres->comp_error(amu_TM_bound_1));

        amu_OS_bound_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound}_(OS)_1",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_OS_bound_1, Njack, jack_file);
        printf("amu_OS(bound)_1 = %g  %g \n", amu_OS_bound_1[Njack - 1], myres->comp_error(amu_OS_bound_1));

        //////////////////////////////////////////////////////////////
        // bounding meff_t
        //////////////////////////////////////////////////////////////
        ibound = 3;

        amu_TM_bound_meff_t_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff_t}_(TM)_1",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_TM_bound_meff_t_1, Njack, jack_file);
        printf("amu_TM(bound_meff_t)_1 = %g  %g \n", amu_TM_bound_meff_t_1[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_1));

        amu_OS_bound_meff_t_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff_t}_(OS)_1",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_OS_bound_meff_t_1, Njack, jack_file);
        printf("amu_OS(bound_meff_t)_1 = %g  %g \n", amu_OS_bound_meff_t_1[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_1));

        //////////////////////////////////////////////////////////////
        // bounding meff
        //////////////////////////////////////////////////////////////
        ibound = 1;

        double* M_VKVK_TM_1 = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack_1, Njack,
            namefile_plateaux, outfile, 0, "M_VKVK_TM_1", M_eff_T, jack_file);

        double* M_VKVK_OS_1 = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack_1, Njack,
            namefile_plateaux, outfile, 1, "M_VKVK_OS_1", M_eff_T, jack_file);

        amu_TM_bound_meff_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff}_(TM)_1",
            resampling, isub, ibound, Mpi, M_VKVK_TM_1, head.L);
        write_jack(amu_TM_bound_meff, Njack, jack_file);
        printf("amu_{bound_meff}_(TM)_1 = %g  %g \n", amu_TM_bound_meff[Njack - 1], myres->comp_error(amu_TM_bound_meff));

        amu_OS_bound_meff_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_{bound_meff}_(OS)_1",
            resampling, isub, ibound, Mpi, M_VKVK_OS_1, head.L);
        write_jack(amu_OS_bound_meff, Njack, jack_file);
        printf("amu_{bound_meff}_(OS)_1 = %g  %g \n", amu_OS_bound_meff[Njack - 1], myres->comp_error(amu_OS_bound_meff));


    }
    else {
        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);

        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);
        zero_corr(zeros, Njack, jack_file);
    }

    double* amu_TM_bound_ave = (double*)malloc(sizeof(double) * Njack);;
    double* amu_OS_bound_ave = (double*)malloc(sizeof(double) * Njack);
    double* amu_TM_bound_meff_t_ave = (double*)malloc(sizeof(double) * Njack);
    double* amu_OS_bound_meff_t_ave = (double*)malloc(sizeof(double) * Njack);
    double* amu_TM_bound_meff_ave = (double*)malloc(sizeof(double) * Njack);
    double* amu_OS_bound_meff_ave = (double*)malloc(sizeof(double) * Njack);
    if (strcmp(option[6], "C80") == 0) {
        weighted_average(amu_TM_bound_ave, amu_TM_bound, amu_TM_bound_1);
        weighted_average(amu_OS_bound_ave, amu_OS_bound, amu_OS_bound_1);
        weighted_average(amu_TM_bound_meff_t_ave, amu_TM_bound_meff_t, amu_TM_bound_meff_t_1);
        weighted_average(amu_OS_bound_meff_t_ave, amu_OS_bound_meff_t, amu_OS_bound_meff_t_1);
        weighted_average(amu_TM_bound_meff_ave, amu_TM_bound_meff, amu_TM_bound_meff_1);
        weighted_average(amu_OS_bound_meff_ave, amu_OS_bound_meff, amu_OS_bound_meff_1);

    }
    else {
        myres->copy(amu_TM_bound_ave, amu_TM_bound);
        myres->copy(amu_OS_bound_ave, amu_OS_bound);
        myres->copy(amu_TM_bound_meff_t_ave, amu_TM_bound_meff_t);
        myres->copy(amu_OS_bound_meff_t_ave, amu_OS_bound_meff_t);
        myres->copy(amu_TM_bound_meff_ave, amu_TM_bound_meff);
        myres->copy(amu_OS_bound_meff_ave, amu_OS_bound_meff);

    }
    write_jack(amu_TM_bound_ave, Njack, jack_file); check_correlatro_counter(22);
    write_jack(amu_OS_bound_ave, Njack, jack_file); check_correlatro_counter(23);
    write_jack(amu_TM_bound_meff_t_ave, Njack, jack_file); check_correlatro_counter(24);
    write_jack(amu_OS_bound_meff_t_ave, Njack, jack_file); check_correlatro_counter(25);
    write_jack(amu_TM_bound_meff_ave, Njack, jack_file); check_correlatro_counter(26);
    write_jack(amu_OS_bound_meff_ave, Njack, jack_file); check_correlatro_counter(27);

    printf("amu_TM(bound)_ave = %g  %g \n", amu_TM_bound_ave[Njack - 1], myres->comp_error(amu_TM_bound_ave));
    printf("amu_OS(bound)_ave = %g  %g \n", amu_OS_bound_ave[Njack - 1], myres->comp_error(amu_OS_bound_ave));
    printf("amu_TM(bound_meff_t)_ave = %g  %g \n", amu_TM_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_ave));
    printf("amu_OS(bound_meff_t)_ave = %g  %g \n", amu_OS_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_ave));
    printf("amu_TM(bound_meff)_ave = %g  %g \n", amu_TM_bound_meff_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_ave));
    printf("amu_OS(bound_meff)_ave = %g  %g \n", amu_OS_bound_meff_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_ave));
    ///// to get the result in quarto
    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_ave[Njack - 1], myres->comp_error(amu_TM_bound_ave));

    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_ave[Njack - 1], myres->comp_error(amu_OS_bound_ave));

    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound_meff_t}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_ave));

    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound_meff_t}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_ave));

    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound_meff}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_meff_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_ave));

    fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
    fprintf(outfile, "\n\n #amu_{bound_meff}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_meff_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_ave));


    //////////////////////////////////////////////////////////////
    // unblind combo
    //////////////////////////////////////////////////////////////
    struct fit_type fit_info;
    fit_info.Njack = Njack;
    fit_info.T = head.T;
    fit_info.corr_id = { 0,1 };

    double* ones = myres->create_fake(1, 0.0000001, 1);
    double* X = unblind_combo(conf_jack, ZA, ZV, fit_info);
    printf("X(%s)= %g   %g \n", argv[3], X[Njack - 1], myres->comp_error(X));
    if (strcmp(option[6], "C80") == 0) {
        double* X1 = unblind_combo(conf_jack_1, ZA, ZV, fit_info);
        printf("X430(%s)= %g   %g \n", argv[3], X1[Njack - 1], myres->comp_error(X1));
        double* X_ave = (double*)malloc(sizeof(double) * Njack);
        weighted_average(X_ave, X, X1);
        printf("X_ave(%s)= %g   %g \n", argv[3], X_ave[Njack - 1], myres->comp_error(X_ave));
        free(X_ave);
        free(X1);
    }
    //////////////////////////////////////////////////////////////
    // FVE GS
    //////////////////////////////////////////////////////////////
    char name_GS[NAMESIZE];
    mysprintf(name_GS, NAMESIZE, "%s/out/%s_GS", option[3], option[6]);
    double* DV;
    double* DV_M;
    double* DV_Mpi_phys;
    double* DV_DM;
    double* DV_C80;
    double* DV_to_C80;
    if (std::filesystem::exists(name_GS)) {
        printf("reading GS from file %s\n", name_GS);
        DV = new double[Njack];
        DV_M = new double[Njack];
        DV_Mpi_phys = new double[Njack];
        DV_DM = new double[Njack];
        DV_C80 = new double[Njack];
        DV_to_C80 = new double[Njack];

        FILE* f = open_file(name_GS, "r+");
        int fi = 00;
        int read_j = -1;
        fi += fscanf(f, "%d\n", &read_j);
        error(read_j != Njack, 1, "main", "Njack in file %g");
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV[j]);
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV_M[j]);
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV_Mpi_phys[j]);
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV_DM[j]);
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV_C80[j]);
        for (int j = 0; j < Njack;j++)    fi += fscanf(f, "%lf\n", &DV_to_C80[j]);
        fclose(f);
        printf("reading GS complete\n");
    }
    else {
        double* Mpi_phys_Mev = myres->create_fake(v_MpiMeV, err_MpiMeV, 1003);
        double* MpiMev = new double[Njack];
        for (int j = 0;j < Njack;j++) {
            MpiMev[j] = Mpi[j] / (a_fm[j] / 197.326963);
        }
        double* jack_Mrho_MeV_exp = fake_sampling(resampling, Mrho_MeV, Mrho_MeV_err, Njack, 1001);
        double* jack_grhopipi = fake_sampling(resampling, grhopipi, grhopipi_err, Njack, 1002);
        DV = compute_DVt_and_integrate(head.L, Njack, MpiMev /* jack_Mpi_MeV_exp */, jack_Mrho_MeV_exp, a_fm, jack_grhopipi, outfile, "DVt", resampling, 24, 1);


        DV_M = compute_DVt_and_integrate(head.L, Njack, MpiMev /* jack_Mpi_MeV_exp */, jack_Mrho_MeV_exp, a_fm, jack_grhopipi, outfile, "DVt_M", resampling, 24, 2);
        DV_Mpi_phys = compute_DVt_and_integrate(head.L, Njack, Mpi_phys_Mev /* jack_Mpi_MeV_exp */, jack_Mrho_MeV_exp, a_fm, jack_grhopipi, outfile, "DVt_Mpiphys", resampling, 24, 2);
        DV_DM = new double[Njack];
        myres->sub(DV_DM, DV_Mpi_phys, DV_M);

        double* a_fm_C80 = myres->create_fake(0.06821, 0.00013, 5);
        int L_C80 = 80;
        DV_C80 = compute_DVt_and_integrate(L_C80, Njack, MpiMev /* jack_Mpi_MeV_exp */, jack_Mrho_MeV_exp, a_fm_C80, jack_grhopipi, outfile, "DVtC80", resampling, 24, 1);

        DV_to_C80 = myres->create_copy(DV);
        myres->sub(DV_to_C80, DV, DV_C80);
        FILE* f = open_file(name_GS, "w+");
        int fi = 00;
        fi += fprintf(f, "%d\n", Njack);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV[j]);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV_M[j]);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV_Mpi_phys[j]);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV_DM[j]);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV_C80[j]);
        for (int j = 0; j < Njack;j++)    fi += fprintf(f, "%.12g\n", DV_to_C80[j]);
        fclose(f);
    }
    printf("GS= %g  %g\n", (10.0 / 9.0) * DV[Njack - 1], (10.0 / 9.0) * myres->comp_error(DV));
    printf("GS_M= %g  %g\n", (10.0 / 9.0) * DV_DM[Njack - 1], (10.0 / 9.0) * myres->comp_error(DV_DM));
    printf("GS_C80= %g  %g\n", (10.0 / 9.0) * DV_C80[Njack - 1], (10.0 / 9.0) * myres->comp_error(DV_C80));
    printf("GS_to_C80= %g  %g\n", (10.0 / 9.0) * DV_to_C80[Njack - 1], (10.0 / 9.0) * myres->comp_error(DV_to_C80));
    write_jack(DV, Njack, jack_file); check_correlatro_counter(28);
    write_jack(DV_DM, Njack, jack_file); check_correlatro_counter(29);
    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < 2; t++) {
        fprintf(outfile, "%d   1   1\t", t);
    }
    fprintf(outfile, "\n\n #DV_DM fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", DV_DM[Njack - 1], error_jackboot(resampling, Njack, DV_DM));

    write_jack(DV_C80, Njack, jack_file); check_correlatro_counter(30);
    write_jack(DV_to_C80, Njack, jack_file); check_correlatro_counter(31);
    fprintf(outfile, " \n\n");
    fprintf(outfile, "#\n");
    for (int t = 1; t < 2; t++) {
        fprintf(outfile, "%d   1   1\t", t);
    }
    fprintf(outfile, "\n\n #DV_to_C80 fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 2, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %15.g\n", DV_to_C80[Njack - 1], error_jackboot(resampling, Njack, DV_to_C80));

    //////////////////////////////////////////////////////////////
    // lattice artefact three sub
    //////////////////////////////////////////////////////////////
    {
        //////////////////////////////////////////////////////////////
        // bounding BMW
        //////////////////////////////////////////////////////////////
        int ibound = 2;

        isub = 2;
        double* amu_TM_bound = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound}_(TM)",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_TM_bound, Njack, jack_file);
        printf("amu_TM(bound) = %g  %g \n", amu_TM_bound[Njack - 1], myres->comp_error(amu_TM_bound));

        isub = 3;
        double* amu_OS_bound = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound}_(OS)",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_OS_bound, Njack, jack_file);
        printf("amu_OS(bound) = %g  %g \n", amu_OS_bound[Njack - 1], myres->comp_error(amu_OS_bound));

        //////////////////////////////////////////////////////////////
        // bounding meff_t
        //////////////////////////////////////////////////////////////
        ibound = 3;
        isub = 2;
        double* amu_TM_bound_meff_t = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff_t}_(TM)",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_TM_bound_meff_t, Njack, jack_file);
        printf("amu_TM(bound_meff_t) = %g  %g \n", amu_TM_bound_meff_t[Njack - 1], myres->comp_error(amu_TM_bound_meff_t));
        isub = 3;
        double* amu_OS_bound_meff_t = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff_t}_(OS)",
            resampling, isub, ibound, Mpi, nullptr, head.L);
        write_jack(amu_OS_bound_meff_t, Njack, jack_file);
        printf("amu_OS(bound_meff_t) = %g  %g \n", amu_OS_bound_meff_t[Njack - 1], myres->comp_error(amu_OS_bound_meff_t));

        //////////////////////////////////////////////////////////////
        // bounding meff
        //////////////////////////////////////////////////////////////
        ibound = 1;

        write_jack(zeros, Njack, jack_file);
        write_jack(zeros, Njack, jack_file);

        isub = 2;
        double* amu_TM_bound_meff = compute_amu_bounding(conf_jack, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff}_(TM)",
            resampling, isub, ibound, Mpi, M_VKVK_TM, head.L);
        write_jack(amu_TM_bound_meff, Njack, jack_file);
        printf("amu_TM(bound_meff) = %g  %g \n", amu_TM_bound_meff[Njack - 1], myres->comp_error(amu_TM_bound_meff));
        isub = 3;
        double* amu_OS_bound_meff = compute_amu_bounding(conf_jack, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff}_(OS)",
            resampling, isub, ibound, Mpi, M_VKVK_OS, head.L);
        write_jack(amu_OS_bound_meff, Njack, jack_file);
        printf("amu_OS(bound_meff) = %g  %g \n", amu_OS_bound_meff[Njack - 1], myres->comp_error(amu_OS_bound_meff));

        double* amu_TM_bound_1;
        double* amu_OS_bound_1;
        double* amu_TM_bound_meff_t_1;
        double* amu_OS_bound_meff_t_1;
        double* amu_TM_bound_meff_1;
        double* amu_OS_bound_meff_1;

        if (strcmp(option[6], "C80") == 0) {
            //////////////////////////////////////////////////////////////
            // bounding BMW
            //////////////////////////////////////////////////////////////
            int ibound = 2;
            isub = 2;
            amu_TM_bound_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound}_(TM)_1",
                resampling, isub, ibound, Mpi, nullptr, head.L);
            write_jack(amu_TM_bound_1, Njack, jack_file);
            printf("amu_TM(bound)_1 = %g  %g \n", amu_TM_bound_1[Njack - 1], myres->comp_error(amu_TM_bound_1));
            isub = 3;
            amu_OS_bound_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound}_(OS)_1",
                resampling, isub, ibound, Mpi, nullptr, head.L);
            write_jack(amu_OS_bound_1, Njack, jack_file);
            printf("amu_OS(bound)_1 = %g  %g \n", amu_OS_bound_1[Njack - 1], myres->comp_error(amu_OS_bound_1));

            //////////////////////////////////////////////////////////////
            // bounding meff_t
            //////////////////////////////////////////////////////////////
            ibound = 3;
            isub = 2;
            amu_TM_bound_meff_t_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff_t}_(TM)_1",
                resampling, isub, ibound, Mpi, nullptr, head.L);
            write_jack(amu_TM_bound_meff_t_1, Njack, jack_file);
            printf("amu_TM(bound_meff_t)_1 = %g  %g \n", amu_TM_bound_meff_t_1[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_1));
            isub = 3;
            amu_OS_bound_meff_t_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff_t}_(OS)_1",
                resampling, isub, ibound, Mpi, nullptr, head.L);
            write_jack(amu_OS_bound_meff_t_1, Njack, jack_file);
            printf("amu_OS(bound_meff_t)_1 = %g  %g \n", amu_OS_bound_meff_t_1[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_1));

            //////////////////////////////////////////////////////////////
            // bounding meff
            //////////////////////////////////////////////////////////////
            ibound = 1;

            double* M_VKVK_TM_1 = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack_1, Njack,
                namefile_plateaux, outfile, 0, "M_VKVK_TM_1", M_eff_T, jack_file);

            double* M_VKVK_OS_1 = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack_1, Njack,
                namefile_plateaux, outfile, 1, "M_VKVK_OS_1", M_eff_T, jack_file);
            isub = 2;
            amu_TM_bound_meff_1 = compute_amu_bounding(conf_jack_1, 0, Njack, ZA, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff}_(TM)_1",
                resampling, isub, ibound, Mpi, M_VKVK_TM_1, head.L);
            write_jack(amu_TM_bound_meff, Njack, jack_file);
            printf("amu_{bound_meff}_(TM)_1 = %g  %g \n", amu_TM_bound_meff[Njack - 1], myres->comp_error(amu_TM_bound_meff));
            isub = 3;
            amu_OS_bound_meff_1 = compute_amu_bounding(conf_jack_1, 1, Njack, ZV, a_fm, q2ud, int_scheme, outfile, "amu_treesub_{bound_meff}_(OS)_1",
                resampling, isub, ibound, Mpi, M_VKVK_OS_1, head.L);
            write_jack(amu_OS_bound_meff, Njack, jack_file);
            printf("amu_{bound_meff}_(OS)_1 = %g  %g \n", amu_OS_bound_meff[Njack - 1], myres->comp_error(amu_OS_bound_meff));


        }
        else {
            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);

            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);
            zero_corr(zeros, Njack, jack_file);
        }

        double* amu_TM_bound_ave = (double*)malloc(sizeof(double) * Njack);;
        double* amu_OS_bound_ave = (double*)malloc(sizeof(double) * Njack);
        double* amu_TM_bound_meff_t_ave = (double*)malloc(sizeof(double) * Njack);
        double* amu_OS_bound_meff_t_ave = (double*)malloc(sizeof(double) * Njack);
        double* amu_TM_bound_meff_ave = (double*)malloc(sizeof(double) * Njack);
        double* amu_OS_bound_meff_ave = (double*)malloc(sizeof(double) * Njack);
        if (strcmp(option[6], "C80") == 0) {
            weighted_average(amu_TM_bound_ave, amu_TM_bound, amu_TM_bound_1);
            weighted_average(amu_OS_bound_ave, amu_OS_bound, amu_OS_bound_1);
            weighted_average(amu_TM_bound_meff_t_ave, amu_TM_bound_meff_t, amu_TM_bound_meff_t_1);
            weighted_average(amu_OS_bound_meff_t_ave, amu_OS_bound_meff_t, amu_OS_bound_meff_t_1);
            weighted_average(amu_TM_bound_meff_ave, amu_TM_bound_meff, amu_TM_bound_meff_1);
            weighted_average(amu_OS_bound_meff_ave, amu_OS_bound_meff, amu_OS_bound_meff_1);

        }
        else {
            myres->copy(amu_TM_bound_ave, amu_TM_bound);
            myres->copy(amu_OS_bound_ave, amu_OS_bound);
            myres->copy(amu_TM_bound_meff_t_ave, amu_TM_bound_meff_t);
            myres->copy(amu_OS_bound_meff_t_ave, amu_OS_bound_meff_t);
            myres->copy(amu_TM_bound_meff_ave, amu_TM_bound_meff);
            myres->copy(amu_OS_bound_meff_ave, amu_OS_bound_meff);

        }
        write_jack(amu_TM_bound_ave, Njack, jack_file); check_correlatro_counter(48);
        write_jack(amu_OS_bound_ave, Njack, jack_file); check_correlatro_counter(49);
        write_jack(amu_TM_bound_meff_t_ave, Njack, jack_file); check_correlatro_counter(50);
        write_jack(amu_OS_bound_meff_t_ave, Njack, jack_file); check_correlatro_counter(51);
        write_jack(amu_TM_bound_meff_ave, Njack, jack_file); check_correlatro_counter(52);
        write_jack(amu_OS_bound_meff_ave, Njack, jack_file); check_correlatro_counter(53);

        printf("amu_treesub_TM(bound)_ave = %g  %g \n", amu_TM_bound_ave[Njack - 1], myres->comp_error(amu_TM_bound_ave));
        printf("amu_treesub_OS(bound)_ave = %g  %g \n", amu_OS_bound_ave[Njack - 1], myres->comp_error(amu_OS_bound_ave));
        printf("amu_treesub_TM(bound_meff_t)_ave = %g  %g \n", amu_TM_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_ave));
        printf("amu_treesub_OS(bound_meff_t)_ave = %g  %g \n", amu_OS_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_ave));
        printf("amu_treesub_TM(bound_meff)_ave = %g  %g \n", amu_TM_bound_meff_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_ave));
        printf("amu_treesub_OS(bound_meff)_ave = %g  %g \n", amu_OS_bound_meff_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_ave));
        ///// to get the result in quarto
        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_ave[Njack - 1], myres->comp_error(amu_TM_bound_ave));

        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_ave[Njack - 1], myres->comp_error(amu_OS_bound_ave));

        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound_meff_t}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_t_ave));

        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound_meff_t}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_meff_t_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_t_ave));

        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound_meff}_(TM)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_TM_bound_meff_ave[Njack - 1], myres->comp_error(amu_TM_bound_meff_ave));

        fprintf(outfile, " \n\n#\n%d   %.15g   %.15g %d   %.15g   %.15g", 0, 0.0, 0.0, 0, 0.0, 0.0);
        fprintf(outfile, "\n\n #amu_treesub_{bound_meff}_(OS)_ave_eigen fit in [%d,%d] chi2=%.5g  %.5g\n", 1, 1, 0.0, 0.0);
        fprintf(outfile, "   %.15g   %15.g\n", amu_OS_bound_meff_ave[Njack - 1], myres->comp_error(amu_OS_bound_meff_ave));
    }
}