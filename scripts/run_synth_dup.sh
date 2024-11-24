#!/bin/bash
output_file=../results/bench_synth_dup_$(date +"%Y%m%d_%I%M%p").csv

cd build
for ((i = 1 ; i < 10 ; i++ ))
do 
    ./bench_synth_bls 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_ips 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_uls600 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_uls1000 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_rs 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_ls21 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_pdq 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_ls20 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_ska 800 dup | tee -a $output_file
    sleep 10
    ./bench_synth_stdsort 800 dup | tee -a $output_file
    sleep 10
done
cd ..
