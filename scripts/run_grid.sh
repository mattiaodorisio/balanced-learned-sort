#!/bin/bash
output_file=../results/bench_grid_$(date +"%Y%m%d_%I%M%p").csv

cd build
for ((i = 1 ; i < 10 ; i++ ))
do 
    ./bench_grid_bls 800 | tee -a $output_file
    sleep 60
    ./bench_grid_uls 800 | tee -a $output_file
    sleep 60
    ./bench_grid_ls21 800 | tee -a $output_file
    sleep 60
done
cd ..
