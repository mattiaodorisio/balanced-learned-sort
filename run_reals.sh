#!/bin/bash
output_file=../results/bench_real_$(date +"%Y%m%d_%I%M%p").csv

cd build
for ((i = 1 ; i < 10 ; i++ ))
do 
    ./bench_real_bls | tee -a $output_file
    sleep 60
    ./bench_real_ips | tee -a $output_file
    sleep 60
    ./bench_real_uls600 | tee -a $output_file
    sleep 60
    ./bench_real_uls1000 | tee -a $output_file
    sleep 60
    ./bench_real_rs | tee -a $output_file
    sleep 60
    ./bench_real_ls21 | tee -a $output_file
    sleep 60
    ./bench_real_pdq | tee -a $output_file
    sleep 60
    ./bench_real_ls20 | tee -a $output_file
    sleep 60
    ./bench_real_ska | tee -a $output_file
    sleep 60
    ./bench_real_stdsort | tee -a $output_file
    sleep 60
done
cd ..
