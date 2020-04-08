#! /bin/sh
#SBATCH --partition=general --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --time=08:00:00
#SBATCH --job-name=stringdb


source ~/env-keras/bin/activate
#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



path=$1
clf=$2
#gba, n2v_knn
fold=$3


python stringExperiment.py $path $clf $fold
