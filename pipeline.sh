#!/bin/sh

#pre-process obo file to remove comments and typedefs
tail go.obo -n+27 | head -n-73 > go-final.obo



#for all species
speciesList=("yeast" "tomato" "arabidopsis" "celegans")
#for one species
#speciesList=("yeast")

for species in ${speciesList[@]};
do

    tail -n+13 data/$species/annotations/goa_$species.gaf > data/$species/annotations/goa.gaf
done


#parse gene annotation (gaf) files to generate a matrix #proteins x #GO terms
#saved in data/$species/annotations
for species in ${speciesList[@]};
do
    python preProcessAnnotations.py $species $ontology $removeIPI
done


annoPref="P_"
thresholdType="median"

#blast and test sets
for species in ${speciesList[@]};
do

  python setupExperiment.py $species $annoPref $thresholdType;
  python preProcessAnnotations.py swissprot;
  python blast_makeSubsets.py "../experiments/"$species"_"$annoPref$thresholdType"/";

  for fold in {0..4};
  do
    makeblastdb -in "../blast/subsets/"$species"training_fold"$fold".fasta" -dbtype prot;
    #blastp -query "../blast/subsets/"$species"test_fold"$fold".fasta" -db "../blast/subsets/"$species"training_fold"$fold".fasta" -out "../blast/"$species"_fold"$fold".txt" -outfmt 7
    #./myblastp.sh $species $fold
    sbatch myblastp.sh $species $fold

  done;
done

exit


#parameters
ontology="P"
removeIPI="0"

thresholdType="median"
thresholdValue=0.5
if [[ $removeIPI == "0" ]]
then
  annoPrefix=$ontology"_"
else
  annoPrefix=$ontology"_noIPI_"
fi


cd code/


python setupExperiment.py $species $annoPrefix $thresholdType $thresholdValue
python blast_makeSubsets.py "../experiments/"$species"_"$ontology$annoPrefix$thresholdType"/"

export PATH=$PATH:/tudelft.net/staff-bulk/ewi/insy/DBL/smakrod/local-blast-install/ncbi-blast-2.5.0+/bin/

for fold in {0..4};
do
  makeblastdb -in "../blast/subsets/"$species"training_fold"$fold".fasta" -dbtype prot;
  sbatch myblastp.sh $species $fold;
  #use this if not in cluster
  #./myblastp.sh $species $fold
done


#parse biogrid files to generate adjacency matrix
for species in ${speciesList[@]};
do

    ./preprocessBiogrid.sh $species
done
