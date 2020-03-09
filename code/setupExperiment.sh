#! /bin/sh
#SBATCH --partition=general --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=50000
#SBATCH --time=08:30:00
#SBATCH --job-name=setup

source ~/env-keras/bin/activate

#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



species=$1
annoPrefix=$2
thresType=$3

python setupExperiment.py $species $annoPrefix $thresType
