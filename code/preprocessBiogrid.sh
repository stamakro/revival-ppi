

path="../data/"$1"/interactions/"
mkdir -p $path"final/biogrid/"

if [ $1 == "yeast" ];
then

	#remove header
	tail -n+36 $path"ppi-biogrid.txt" > $path"ppi-biogrid-noheader.txt";
	#remove interactions with proteins from other species
	awk -F "\t" '$10 == "559292" && $11=="559292" {print $0}' $path"ppi-biogrid-noheader.txt" > $path"ppi-biogrid-noheader-withinspecies.txt";

	#keep physical interactions only
	grep -v -f biogrid-genetic-interactions-codes.txt $path"ppi-biogrid-noheader-withinspecies.txt" > $path"ppi-biogrid-physical.txt";

	python makeNetwork.py $1;

	python createAdjacencyMatrix.py $1;

elif [ $1 == "arabidopsis" ];
then
	#remove header
	tail -n+36 $path"ppi-biogrid.txt" > $path"ppi-biogrid-noheader.txt";

	#remove interactions with proteins from other species
	awk -F "\t" '$10 == "3702" && $11=="3702" {print $0}' $path"ppi-biogrid-noheader.txt" > $path"ppi-biogrid-noheader-withinspecies.txt";

	#keep physical interactions only
	grep -v -f biogrid-genetic-interactions-codes.txt $path"ppi-biogrid-noheader-withinspecies.txt" > $path"ppi-biogrid-physical.txt";

	python makeNetwork.py $1;

	python createAdjacencyMatrix.py $1;



elif [ $1 == "tomato" ];
then
	#remove header
	tail -n+36 $path"ppi-biogrid.txt" > $path"ppi-biogrid-noheader.txt";
	#remove interactions with proteins from other species
	awk -F "\t" '$10 == "4081" && $11=="4081" {print $0}' $path"ppi-biogrid-noheader.txt" > $path"ppi-biogrid-noheader-withinspecies.txt";
	#keep physical interactions only
	grep -v -f biogrid-genetic-interactions-codes.txt $path"ppi-biogrid-noheader-withinspecies.txt" > $path"ppi-biogrid-physical.txt";

	python makeNetwork.py $1;

	python createAdjacencyMatrix.py $1;

elif [ $1 == "celegans" ];
then
	#remove header
	tail -n+36 $path"ppi-biogrid.txt" > $path"ppi-biogrid-noheader.txt";
	#remove interactions with proteins from other species
	awk -F "\t" '$10 == "6239" && $11=="6239" {print $0}' $path"ppi-biogrid-noheader.txt" > $path"ppi-biogrid-noheader-withinspecies.txt";
	#keep physical interactions only
	grep -v -f biogrid-genetic-interactions-codes.txt $path"ppi-biogrid-noheader-withinspecies.txt" > $path"ppi-biogrid-physical.txt";

	python makeNetwork.py $1;

	python createAdjacencyMatrix.py $1;




fi
