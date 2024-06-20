#!/bin/bash
#$ -cwd
#$ -N ASSR_CochTune
#$ -q cpu.q
#$ -pe smp 60
<<<<<<< HEAD
#$ -l h_vmem=256G
#$ -l h_rt=1:30:00
#$ -o /ddn/yubelyrn/A1model_sm/data/singleSim.out
#$ -e /ddn/yubelyrn/A1model_sm/data/singleSim.err
=======
>>>>>>> upstream/main

source ~/.bashrc
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi init.py