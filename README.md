# Kython

Kmer Computation tools python library

This package was developped for a university project at Sorbonne université for the use of Kmer signatures in phylogenetic studies.


### Installation:

```bash

pip3 install kython

```

### Functions available:

* Module pp:

        * DownloadSequences(Bacteria list file, Archaea list file, output path)-Download sequences From ncbi and returns a dictionnary of file paths
        * ParseSequences(path)--------------------------------------------------Retrieve the file path dictionnary from an existing folder
        * KmerSignature(path,kmer size,normalization) --------------------------Compute Kmer signature from a genomic .fna file path (use the dictionnary)
        * DistanceMatrix(dictGeneral:dict, kmer:int,phylum:str=None)------------Compute distance matrices from a dictionnary (phylum can be specified, default is all phylums combined)
        * NeighbourJoining(distance Matrix)-------------------------------------Compute Phylogenetic tree in Newick format from distance Matrix (using Neigbour Joining algorithm)
        * SequenceHomogeneity(path, kmer size, Genomic fragment size)-----------Compute Homgeneiy of sub sequences in the main sequence (using chisquare) and return list of intervals and p-values
        
* Module pl:
        * KmerSignature(signature,organism name)--------------------------------Plots a signature in radial coordinates
        * SequenceHomogeneity(listPvalues,listPositions,pValueAccepted,organism name, kmer size)-Plots the sequence homogeneity
        * GraphTree(tree in newick format,DictGeneral,outputpath)
