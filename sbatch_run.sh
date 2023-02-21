#!/bin/bash

num_proc=2000

sbatch --time 167:00:00 --ntasks $num_proc --qos standby --requeue --mem-per-cpu=32000M --mail-user jastern33@gmail.com --mail-type BEGIN --mail-type END --mail-type FAIL fragment_screening_wrapper.sh $num_proc
# sbatch --time 00:15:00 --ntasks $num_proc --mem-per-cpu=32000M --mail-user jastern33@gmail.com --mail-type BEGIN --mail-type END --mail-type FAIL fragment_screening_wrapper.sh $num_proc
