#! /bin/sh
#SBATCH --partition=general --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000
#SBATCH --time=06:00:00
#SBATCH --job-name=n2v


#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

inputfile=$1
nwalks=$2
p=$3
q=$4
length=$5

python node2vec/src/main.py --input $inputfile --walk-length $length --workers 2 --p $p --q $q --num-walks $nwalks;
