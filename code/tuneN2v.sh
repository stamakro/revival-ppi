
species="yeast"
inputfile="../data/"$species"/networks/tmp129.network"

for nwalks in 5 10 15 20 30
do
	for p in 0.1 0.5 1.0 1.5 3.0 
	do
		for q in 0.1 0.5 1.0 1.5 3.0 
		do
			for length in 20 40 60 80 100 160
			do
			
				if ! [ -f "node2vec/emb/string/"$species"_"$nwalks"_"$p"_"$q"_"$length"_500.emb" ]
				then
					echo $inputfile $nwalks $p $q $length
					sbatch tuneN2vWorker.sh $inputfile $nwalks $p $q $length
					#./tuneN2vWorker.sh $inputfile $nwalks $p $q $length
					#exit
				fi			
				
			done
		done
	done
done

