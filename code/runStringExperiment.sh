#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --time=02:00:00
#SBATCH --job-name=stringdb


source ~/env-keras/bin/activate
#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



path=$1
clf="gba"
fold=$2

python stringExperiment.py $path $clf $fold
