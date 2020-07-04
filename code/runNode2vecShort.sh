#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=02:25:00
#SBATCH --job-name=n2v


#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



species=$1
inputfile="../data/"$1"/networks/tmp"$2".network"
nwalks=10
p=0.5
q=1.0
length=20
dim=500


mkdir -p ../data/"$1"/networks/"$nwalks"_"$p"_"$q"_"$length"_"$dim"/

python node2vec/src/main.py --input $inputfile --walk-length $length --workers 2 --p $p --q $q --num-walks $nwalks --dimensions $dim --output "../data/"$1"/networks/"$nwalks"_"$p"_"$q"_"$length"_"$dim"/tmp"$2".emb" ;	


#if [ ! -f "../data/"$1"/networks/tmp"$2".emb" ];
#then
#	python node2vec/src/main.py --input "../data/"$1"/networks/tmp"$2".network" --output "../data/"$1"/networks/tmp"$2".emb" 
	
#fi

