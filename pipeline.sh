#!/bin/sh

#pre-process obo file to remove comments and typedefs
if ! [ -f ../data/go/go-final.obo  ]
then
	tail go.obo -n+27 | head -n-73 > go-final.obo;
fi

species="tomato"
ontology="P"

thresholdType="median"
removeIPI="0"

thresholdValue=0.5
if [[ $removeIPI == "0" ]]
then
  annoPrefix=$ontology"_"
else
  annoPrefix=$ontology"_noIPI_"
fi



tail -n+13 data/$species/annotations/goa_$species.gaf > data/$species/annotations/goa.gaf

cd code/

#parse gene annotation (gaf) files to generate a matrix #proteins x #GO terms
#saved in data/$species/annotations
python preProcessAnnotations.py $species $ontology $removeIPI

./preprocessBiogrid.sh $species
./preprocessString.sh $species


#blast and test sets
python setupExperiment.py $species $annoPrefix $thresholdType;
python preProcessAnnotations.py swissprot;
python blast_makeSubsets.py "../experiments/"$species"_"$annoPrefix$thresholdType"/";

export PATH=$PATH:/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/local-blast-install/ncbi-blast-2.5.0+/bin/

for fold in {0..4};
do
    makeblastdb -in "../blast/subsets/"$species"training_fold"$fold".fasta" -dbtype prot;
    #blastp -query "../blast/subsets/"$species"test_fold"$fold".fasta" -db "../blast/subsets/"$species"training_fold"$fold".fasta" -out "../blast/"$species"_fold"$fold".txt" -outfmt 7
    #./myblastp.sh $species $fold
    sbatch myblastp.sh $species $fold

done;

python blast_predict.py "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/"

sbatch runNode2vec.sh $species
#./runNode2vec.sh $species

for fold in {0..4};
do
	sbatch runStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" "gba" $fold;
	#./runStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" "gba" $fold;

done


for fold in {0..4};
do
	sbatch runStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" "n2v_knn" $fold;
	#./runStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" "n2v_knn" $fold;

done

for fold in {0..4};
do
	for hops in 1 2;
	do
	sbatch evalStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" $hops $fold
	#./evalStringExperiment.sh "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/" $hops $fold
	done
done

python naive.py "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/"
python plotPerformance.py "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/"


exit




