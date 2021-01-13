import subprocess
import time
import os
import re
import itertools
import pandas as pd
from Bio import SeqIO
from scipy.stats import chi2_contingency
from scipy.spatial import distance

""" First Function Downloading the genomes """

def Cleaning_Folder(path):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.gz'):
                subprocess.run("gunzip "+root+"/"+file,shell=True)
            if file=="MD5SUMS":
                subprocess.run("rm "+root+"/"+file,shell=True)
                
                
                        
#Returns a dictionnary: Key: Phylum -> Value Other dictionnary | Organism name -> Path to Proteome and Genome

def Parsing_Sequences(path):
    #Dictionnary of all the adresses
    dictGen={}
    #Getting the species names via Regular expression
    regex="[0-9] (.+)(chromosome|,)"
    
    listPhylums=os.listdir(path)
    listFiles=[[] for _ in range(len(listPhylums))]
    for i in range(len(listPhylums)):
        listfaa=[]
        listfna=[]
        listNames=[]
        for root,dirs,files in os.walk(path+'/'+listPhylums[i]):
            for file in files:
                if file.endswith('.faa'):
                    listfaa.append(root+'/'+file)
                elif file.endswith('.fna'):
                    listfna.append(root+'/'+file)
                    lineSpecie=open(root+'/'+file).readlines()[0]
                    match=re.search(regex,lineSpecie).group(1)
                    if match.split(' ')[-1]=='chromosome':
                        match=' '.join(match.split(' ')[:-1])
                    listNames.append(match)    
        dictGen[listPhylums[i]]=dict(zip(listNames,zip(listfaa,listfna)))          
        
    return dictGen
    
    
def List_Missing_Organisms(bacteriaPath, archaeaPath,filePath):

        dictGeneral=Parsing_Sequences(filePath)
        
        print(dictGeneral)

        listBacterias=[i[:-1] for i in open(bacteriaPath,'r').readlines()[1:]]
        listArchaeas=[i[:-1] for i in open(archaeaPath,'r').readlines()[1:]]
        
        print('Following Bacterias Genomes were not downloaded:','\n')
        for i in listBacterias:
                if i not in dictGeneral['bacteria'].keys():
                        print(i)
        
        print('Following Archaeas Genomes were not downloaded:','\n')    
        for i in listArchaeas:
                if i not in dictGeneral['archaea'].keys():
                        print(i)	
                        

def DownloadingSequences(bacteriaPath,archaeaPath,outputPath):

        timeInit=time.time()

        print("Downloading Files...",'\n')

        print("Downloading Bacteria Files")
        subprocess.run("ncbi-genome-download --genera \'"+bacteriaPath+"\' bacteria -F protein-fasta,fasta -o "+outputPath,shell=True)

        print("Downloading Archaea Files")
        subprocess.run("ncbi-genome-download --genera \'"+archaeaPath+"\' archaea -F protein-fasta,fasta -o "+outputPath,shell=True)

        Cleaning_Folder(outputPath)
        List_Missing_Organisms(bacteriaPath,archaeaPath,outputPath)	

        print('Elapsed time:',time.time()-timeInit)
        
        return True

""" Parsing all genomes """

def ParsingSequences(path):
    #Dictionnary of all the adresses
    dictGen={}
    #Getting the species names via Regular expression
    regex="[0-9] (.+)(chromosome|,)"
    
    listPhylums=os.listdir(path)
    listFiles=[[] for _ in range(len(listPhylums))]
    for i in range(len(listPhylums)):
        listfaa=[]
        listfna=[]
        listNames=[]
        for root,dirs,files in os.walk(path+'/'+listPhylums[i]):
            for file in files:
                if file.endswith('.faa'):
                    listfaa.append(root+'/'+file)
                elif file.endswith('.fna'):
                    listfna.append(root+'/'+file)
                    lineSpecie=open(root+'/'+file).readlines()[0]
                    match=re.search(regex,lineSpecie).group(1)
                    if match.split(' ')[-1]=='chromosome':
                        match=' '.join(match.split(' ')[:-1])
                    listNames.append(match)    
        dictGen[listPhylums[i]]=dict(zip(listNames,zip(listfaa,listfna)))          
        
    return dictGen

""" Retrieve signature """

def Read_Sequence(path):
        seqs=SeqIO.parse(path,'fasta')
        seqs=[str(seq.seq) for seq in seqs]
        return ''.join(seqs)

def Count_Cuts(listOfSequences,normalized=False):
        #Creating the dictionnary
        possibilities=list(map(''.join,list(itertools.product('ACGT', repeat=len(listOfSequences[0])))))
        counts=[0 for i in range(len(possibilities))]
        dicoCuts=dict(zip(possibilities,counts))
        #Counting sequences
        for sequence in listOfSequences:
                try:
                        dicoCuts[sequence]+=1
                except:
                        None
        #Conversion to df
        df=pd.DataFrame([dicoCuts])

        if normalized==False:        
                return df
        else:
                return df/np.sum(df.values)
                
def KmerSignature(path,kmer,normalized):
	sequence = Read_Sequence(path)
	seqCut = [sequence[i:i+kmer] for i in range(len(sequence)-(kmer-1)) ]
	dicKmer = Count_Cuts(seqCut,normalized)
	return dicKmer
	
""" Distance matrix """

def DistanceMatrix(path,kmer):
        start = time.time()
        pathGenomes=path
        dictGeneral=ParsingSequences(pathGenomes)
        matrice=[]
        for i in dictGeneral[filum]:
                pathTest=(dictGeneral[filum][i])[1]
                dicSeq=KmerSignature(pathTest,kmer,normalized=True)
                matrice.append(dicSeq)
        matrice_Distance=np.zeros((len(matrice),len(matrice)))
        for i in range(len(matrice)):
                for j in range(i,len(matrice)):
                        if i!=j:
                                a=matrice[i].values[0]
                                b=matrice[j].values[0]
                                dst = distance.euclidean(a, b)
                                matrice_Distance[i][j]=dst
                                matrice_Distance[j][i]=dst

        matrice_distance_df=pd.DataFrame(data=matrice_Distance,columns=list(dictGeneral[filum].keys()),index=list(dictGeneral[filum].keys()))
        return matrice_distance_df



