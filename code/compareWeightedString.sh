#! /bin/sh
#SBATCH --partition=general --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=40000
#SBATCH --time=20:00:00
#SBATCH --job-name=stringdb


source ~/env-keras/bin/activate
#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


python compareWeightedString.py $1
