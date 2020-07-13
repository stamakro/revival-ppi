#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=01:30:00
#SBATCH --job-name=n2v

source ~/env-keras/bin/activate

#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

python node2vecClassify.py $@
