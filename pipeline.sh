#!/bin/sh

#for all species
speciesList=("yeast" "tomato" "arabidopsis")
#for one species
#speciesList=("yeast")


#parameters
ontology="P"
removeIPI="0"


cd code/

#parse biogrid files to generate adjacency matrix
for species in ${speciesList[@]};
do
    
    ./preprocessBiogrid.sh $species
done


exit

#parse gene annotation (gaf) files to generate a matrix #proteins x #GO terms
#saved in data/$species/annotations
for species in ${speciesList[@]};
do
    python preProcessAnnotations.py $species $ontology $removeIPI
done
