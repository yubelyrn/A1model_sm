#!/bin/bash
#$ -cwd
#$ -N ASSR_CochTune
#$ -q cpu.q
#$ -pe smp 60
#$ -l h_vmem=256G
#$ -l h_rt=1:30:00
#$ -o /ddn/yubelyrn/A1model_forked/data/singleSim.out
#$ -e /ddn/yubelyrn/A1model_forked/data/singleSim.err

source ~/.bashrc
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi init.py