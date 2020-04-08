#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=02:00:00
#SBATCH --job-name=blastp

source ~/env-keras/bin/activate
export PATH=$PATH:/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/local-blast-install/ncbi-blast-2.5.0+/bin/


species=$1
fold=$2


blastp -query "../blast/subsets/"$species"test_fold"$fold".fasta" -db "../blast/subsets/"$species"training_fold"$fold".fasta" -out "../blast/"$species"_fold"$fold".txt" -outfmt 7

