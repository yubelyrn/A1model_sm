#!/bin/bash
#$ -cwd
#$ -N control_test
#$ -q cpu.q
#$ -pe smp 60
#$ -l h_vmem=256G
#$ -l h_rt=3:00:00
#$ -o /ddn/smcelroy97/A1model_sm/control_test.out
#$ -e /ddn/smcelroy97/A1model_sm/control_test.err

source ~/.bashrc
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi init.py