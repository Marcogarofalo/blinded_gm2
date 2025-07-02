./blinded_gm2 -p ../data/ C80 -bin 50 jack  0.00060 || exit 1
./blinded_gm2 -p ../data/ B64 -bin 50 jack 0.00072 || exit 1
./blinded_gm2 -p ../data/ D96 -bin 50 jack 0.00054 || exit 1

./blinded_gm2 -p ../data/ D96_122confs -bin 50 jack 0.00054 || exit 1
./blinded_gm2 -p ../data/ D96_184confs -bin 50 jack 0.00054 || exit 1
./blinded_gm2 -p ../data/ D96_256confs -bin 50 jack 0.00054 || exit 1
./blinded_gm2 -p ../data/ D96_444confs -bin 50 jack 0.00054 || exit 1
./fit_all_blind jack ../data/jackknife/ ../data/fit_all || exit 1


./rew_gm2 -p ../data_tov/ B64 -bin 50 jack 0.00072 reweight_strange_OS_B64_LMA.dat
./rew_gm2 -p ../data_tov/ B64 -bin 50 jack 0.00072 reweight_light_OS_B64_LMA.dat