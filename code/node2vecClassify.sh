species="yeast"
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
					file="../experiments/n2v/"$species"/string/"$fold"_"$nwalks"_"$p"_"$q"_"$length".txt";
					if ! [ -f $file ]
					then
			
						sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length 1;
						
						#echo $fold $nwalks $p $q $length;
						#exit
					else
						n=$(wc -l $file | cut -d " " -f 1);
						if [ $n -lt 70 ]
						then
							sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length 1;
						fi					
					fi			
					continue

					file="../experiments/n2v/"$species"/"$fold"_"$nwalks"_"$p"_"$q"_"$length".txt";
					if ! [ -f $file ]
					then
			
						sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length 0;
						
						#echo $fold $nwalks $p $q $length;
						#exit
					else
						n=$(wc -l $file | cut -d " " -f 1);
						if [ $n -lt 70 ]
						then
							sbatch node2vecClassifyWorker.sh "../experiments/"$species"_P_median/" $fold $nwalks $p $q $length 0;
						fi					
					fi			


				
				done
			done
		done
	done
done
