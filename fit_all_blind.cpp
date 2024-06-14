#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"

// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "resampling_new.hpp"
#include "global.hpp"
#include "do_analysis_charm.hpp"
#include "functions_blinded_gm2.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <map>


enum enum_ensembles_blind_gm2 {
    B64 = 0,
    C80 = 1
};
// generic_header read_header(FILE* stream) {
//     generic_header header;
//     int ir = 0;
//     ir += fread(&header.T, sizeof(int), 1, stream);
//     ir += fread(&header.L, sizeof(int), 1, stream);
//     int s;
//     ir += fread(&s, sizeof(int), 1, stream);
//     header.mus = std::vector<double>(s);
//     for (int i = 0; i < s; i++) {
//         ir += fread(&header.mus[i], sizeof(double), 1, stream);
//     }
//     printf("%d \n",s);
//     exit(1);
//     ir += fread(&s, sizeof(int), 1, stream);
//     header.thetas = std::vector<double>(s);
//     for (int i = 0; i < s; i++) {
//         ir += fread(&header.thetas[i], sizeof(double), 1, stream);
//     }

//     ir += fread(&header.Njack, sizeof(int), 1, stream);
//     header.struct_size = ftell(stream);
//     return header;


// }


double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;

}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;

    //read_single_Njack_Nobs(stream, header.header_size, Njack, Nobs);
    // data_single dj(Nobs,Njack);
    data_single dj;
    dj.header.read_header_jack(stream);
    dj.header.print_header();
    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;

}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    //jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    int count = 0;
    for (std::string s : files) {
        FILE* f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;

}


int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];
    char namefit[NAMESIZE];


    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_B64", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_C80", argv[2], argv[1]);
    files.emplace_back(namefile);

    std::vector<int> myen(files.size());
    for (int e = 0; e < files.size(); e++) {
        myen[e] = e;
        printf("file:  %s\n", files[e].c_str());
    }


    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.create_generalised_resampling();

    std::vector<int> myen_full(jackall.ens);
    for (int e = 0; e < jackall.ens; e++) {
        myen_full[e] = e;
        printf("file:  %s\n", files[e].c_str());
        printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        for (double mu : jackall.en[e].header.mus)
            printf("%g   ", mu);
        printf("\n");
    }
    int Njack = jackall.en[0].Njack;
    std::vector<int> myen_charm = { 0, 2, 3, 4, 5, 6 };

    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }

    double* jack_Mpi_MeV_exp = fake_sampling(argv[1], Mpi_MeV, Mpi_MeV_err, Njack, 1003);

    fit_type fit_info;
    fit_info.N = 2;
    fit_info.Npar = 3;
    fit_info.Nvar = 2;
    fit_info.Njack = Njack;
    fit_info.myen = { B64, C80 };
    fit_info.init_Nxen_from_N_myen();
    fit_info.malloc_x();
    int count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = jackall.en[e].jack[0][j] * jackall.en[e].jack[0][j];  // a^2
                fit_info.x[1][count][j] = jackall.en[e].jack[3][j];  // Mpi
            }
            count++;
        }
    }
    fit_info.corr_id = { 18,19 };
    fit_info.function = rhs_two_line;
    fit_info.linear_fit = true;
    fit_info.covariancey = true;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_identity);
    int ie = 0, ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_bound_a2");
    fit_result amu_linear = fit_all_data(argv, jackall, lhs_identity, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    std::vector<double> xcont = {};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_linear, amu_linear, 0, fit_info.myen.size() - 1, 0.0005, xcont);

}
