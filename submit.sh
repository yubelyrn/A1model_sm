#!/bin/bash
#$ -cwd
#$ -N ASSR_CochTune
#$ -q cpu.q
#$ -pe smp 30
#$ -l h_vmem=256G
#$ -l h_rt=2:00:00
#$ -o /ddn/smcelroy97/A1model_sm/data/singleSim2.out
#$ -e /ddn/smcelroy97/A1model_sm/data/singleSim2.err

source ~/.bashrc
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi init.py