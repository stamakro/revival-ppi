#! /bin/sh
#SBATCH --partition=general --qos=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500
#SBATCH --time=04:00:00
#SBATCH --job-name=n2v

species="tomato"
for fold in {0..4};
do
	for nwalks in 5 10 15 20 30;
	do
		for p in 0.1 0.5 1.0 1.5 3.0;
		do
			for q in 0.1 0.5 1.0 1.5 3.0; 
			do
				for length in 20 40 60 80 100 160;
				do
					file="../experiments/n2v/"$species"/withUnannotated/"$fold"_"$nwalks"_"$p"_"$q"_"$length".txt";
					if ! [ -f $file ]
					then
			
						#sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/withUnannotated/";
						./node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/withUnannotated/";
						
						#echo $fold $nwalks $p $q $length;
						#exit
					else
						n=$(wc -l $file | cut -d " " -f 1);
						if [ $n -lt 70 ]
						then
							#sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/withUnannotated/";
							./node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/withUnannotated/";
						fi					
					fi			
					continue

					file="../experiments/n2v/"$species"/string/"$fold"_"$nwalks"_"$p"_"$q"_"$length".txt";
					if ! [ -f $file ]
					then
			
						sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/string/";
						
						#echo $fold $nwalks $p $q $length;
						#exit
					else
						n=$(wc -l $file | cut -d " " -f 1);
						if [ $n -lt 70 ]
						then
							sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length "/string/";
						fi					
					fi			
					continue

					file="../experiments/n2v/"$species"/"$fold"_"$nwalks"_"$p"_"$q"_"$length".txt";
					if ! [ -f $file ]
					then
			
						sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length;
						
						#echo $fold $nwalks $p $q $length;
						#exit
					else
						n=$(wc -l $file | cut -d " " -f 1);
						if [ $n -lt 70 ]
						then
							sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length;
						fi					
					fi			


				
				done
			done
		done
	done
done
