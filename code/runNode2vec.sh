#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=00:30:00
#SBATCH --job-name=n2v


#so that numpy does not take up all the threads
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1



species=$1

numberNets=$(ls ../data/$1/networks/*.network | wc -l)

echo $numberNets

for (( i=0; i < $numberNets; i++ )); 
do
	echo $i;
	if [ ! -f "../data/"$1"/networks/tmp"$i".emb" ];
	then
		python node2vec/src/main.py --input "../data/"$1"/networks/tmp"$i".network" --output "../data/"$1"/networks/tmp"$i".emb" 
	fi

done
