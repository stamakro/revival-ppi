# PPI Networks with missing edges
This repository contains code for reproducing the article *Sequence-based protein-protein interaction predictions did not outperform protein associations from STRING at functional annotation*

## Summary
1. Download the data
2. Install dependencies 
3. Execute script `pipeline.sh`


## 1. Data
All data used in these experiments for *Sacharomyces cerevisiae*, *Escherichia coli* *Arabidopsis thaliana* and *Solanum lycopersicum* are from the public domain and can be downloaded from the corresponding websites, as detailed below. Data needs to be save at the appropriate directory. The species names are "yeast", "ecoli", "arabidopsis" and "tomato"

### PPI data from BIOGRID
We used version 3.5.181
https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.181/
Save at data/$speciesName$/interactions/ppi-biogrid.txt

### PPI data from STRING
We used version 11.0
https://string-db.org/cgi/download.pl
Save at data/$speciesName$/interactions/ppi-string.txt

### Protein sequences
Protein sequences were downloaded from SwissProt version 2019_11:
ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2019_11/knowledgebase/
Save at data/$speciesName$/sequences/sequences.fasta

For tomato, sequences from the solanum genomics network (https://solgenomics.net/) were also used, available at data/tomato/sequences/sequences_tomatoIDs.fasta

### GO annotations
We used the Gene Ontology released on January 1st 2020 (https://zenodo.org/record/3597529) which is provided at data/go/go.obo

GO annotations of *S. cerevisae* (version 59), *A. thaliana* (version 170) and SwissProt proteins were downloaded from the Gene Ontology Annotation website ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/old/
For tomato, 

## 2. Dependencies
For most analyses, we used Python 3.6 and the packages
+ numpy 1.16.2
+ scipy 1.2.1
+ scikit-learn 0.20.3

For running PIPR, clone https://github.com/muhaochen/seq_ppi
The Keras package is required. We used version  2.2.4

For running node2vec, clone https://github.com/aditya-grover/node2vec
This requires Python 2.7 with networkx 2.2 and numpy 1.14.3

## 3. Execution
After completing all previous steps, execute the following commands:
`./pipeline.sh yeast`
`./pipeline.sh ecoli`
`./pipeline.sh arabidopsis`
`./pipeline.sh tomato`

Note that several steps (running BLAST, training of PIPR, tuning the hyperparameters of node2vec or evaluating all possible combinations of STRING networks) can take several hours to days unless run massively in parallel. All the experiments described in the article were conducted on a compute cluster multiple CPU's and GPU's. 
Duration varies a lot depending on the species. Tomato is the smallest, so it can be run in a few hours. Repeating all experiments for arabidopsis without parallelization will take more than one week.

Run `drawFinalFigures.py` to obtain Figure 2 (with values hard-coded, no experiments needed).
Use `analyzePerformance.py` to measure the effect of individual STRING data sources on the total performances.
Use `plotPerformance.py` to plot the performance and coverage of STRING data sources and combinations thereof.
Use `comparePerProtein.py` to plot performance as a function of node degree (as in Figure 3). See inside this script for possible combinations.
