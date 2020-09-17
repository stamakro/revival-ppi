
species=$1

mkdir -p "../experiments/n2v/"$species"/withUnannotated/"

#create new network
python makeNetwork.py $species "P_" 0
python createAdjacencyMatrix.py $species "P_" 0

python makeNetworkString.py $species "P_" 0
#python createAdjacencyMatrixString.py $species "P_" 0


#convert network to tab delimited node2vec format
python convert2node2vecFormat.py $species

./tuneNode2vec.sh




