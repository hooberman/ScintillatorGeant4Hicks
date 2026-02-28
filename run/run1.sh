cd ../build
make
cd -
rm -f output/sipm_hits.txt
../build/sim_bc408_1

#../build/sim_bc408_1 -n 1 -l "test2Rounded"
