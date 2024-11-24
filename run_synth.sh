#!/bin/bash
output_file=../results/bench_synth_$(date +"%Y%m%d_%I%M%p").csv

cd build
for ((i = 1 ; i < 10 ; i++ ))
do 
    ./bench_synth_bls 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_ips 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_uls 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_rs 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_ls21 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_pdq 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_ls20 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_ska 800 nodup | tee -a $output_file
    sleep 60
    ./bench_synth_stdsort 800 nodup | tee -a $output_file
    sleep 60
done
cd ..
