

path="../data/"$1"/interactions/";
mkdir -p $path"final/string/";
mkdir -p $path"features/";

python makeNetworkString.py $1;
python createAdjacencyMatrixString.py $1;

