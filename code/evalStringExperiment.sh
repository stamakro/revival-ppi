#! /bin/sh
#SBATCH --partition=general --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=02-00:00:00
#SBATCH --job-name=evalstring
#SBATCH --mail-user=stavrosmakrodi
#SBATCH --mail-type=ALL

source ~/env-keras/bin/activate
#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



path=$1
hops=$2
fold=$3
combo=$4

python evaluateGBA.py $path $hops $fold $combo
