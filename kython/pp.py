import subprocess
import time
import os
import re
import itertools
import pandas as pd
import numpy as np
import ncbi_genome_download as ngd
from Bio import SeqIO
from scipy.stats import chisquare
from scipy.spatial import distance
from typing import Tuple



""" Functions written to be used directly by the user:

- DownloadSequences
- KmerSignature
- DistanceMatrix
- NeighbourJoining




"""

""" First Function Downloading the genomes """

def Cleaning_Folder(path:str)->None:
        """Cleans the folder of all .gz file archives and download recap MD5SUMS"""
        for root, dirs, files in os.walk(path):
                for file in files:
                        if file.endswith('.gz'):
                                subprocess.run("gunzip "+root+"/"+file+" -f",shell=True)
                        if file=="MD5SUMS":
                                subprocess.run("rm "+root+"/"+file,shell=True)
                
                
                        
#Returns a dictionnary: Key: Phylum -> Value Other dictionnary | Organism name -> Path to Proteome and Genome

def ParseSequences(path:str) -> dict:
    """ Retrieve a dictionnary from a folder containing phylums folder and then organisms folder
    ex: 
    *refseq
        *bacterias
                *organism1
                *organism2
        *archaeas
                *organism1
                *organism2 """
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
                    
        dictGen[listPhylums[i]]=dict(zip(listNames,listfna))         
        
    return dictGen
    
    
def List_Missing_Organisms(bacteriaPath, archaeaPath,filePath):

        dictGeneral=ParseSequences(filePath)
        
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
                        
        return dictGeneral	
                        

def DownloadSequences(bacteriaPath: str,archaeaPath: str,outputPath: str) -> dict:
        """Download Bacterias Genomes and Archeas Genoms from a list of species and return a dictionnary who links organisms names to genomes files paths"""

        bacterias=open(bacteriaPath,'r').readlines()
        archaeas=open(archaeaPath,'r').readlines()
        
        bacterias=[bacteria[:-1] if bacteria[-1]=='\n' else bacteria for bacteria in bacterias if bacteria!='\n']
        archaeas= [archaea[:-1]  if archaea[-1]=='\n'  else archaea  for archaea in archaeas if archaea!='\n']
        
        timeInit=time.time()

        print("Downloading Files...",'\n')

        print("Downloading Bacteria Files")
        
        for bacteria in bacterias:
                print("Downloading:",bacteria)
                try:
                        ngd.download(section='refseq', 
                                     file_formats='fasta',
                                     genera=bacteria, 
                                     groups='bacteria',
                                     output=outputPath)
                except:
                        print(bacteria+" was not found on NCBI'")

        print("Downloading Archaea Files")
        
        for archaea in archaeas:
                print("Downloading:",archaea)
                try:
                        ngd.download(section='refseq', 
                                     file_formats='fasta',
                                     genera=archaea,
                                     groups='archaea',
                                     output=outputPath)
                except:
                        print(archaea+" was not found on NCBI'")

        Cleaning_Folder(outputPath)
        dictGeneral=List_Missing_Organisms(bacteriaPath,archaeaPath,outputPath+'/refseq/')	

        print('Elapsed time:',time.time()-timeInit)
        
        return dictGeneral


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
                
def KmerSignature(path: str,kmer: int,normalized: bool) -> pd.DataFrame:
        """Computes Kmer Signature from a Genomic .fna file (Use ParseSequences function to retrieve all the adresses from a refseq folder if they are already downloaded, else use DownloadSequences first)"""
        sequence = Read_Sequence(path)
        seqCut = [sequence[i:i+kmer] for i in range(len(sequence)-(kmer-1)) ]
        dicKmer = Count_Cuts(seqCut,normalized)
        return dicKmer





def DistanceMatrix(dictGeneral:dict, kmer:int,phylum:str=None) -> pd.DataFrame:
        """Computes Distance matrix from a dictionnary of file paths (see ParseSequences), if phylum is specified, compute distance matrix only for specified phylum"""
        start = time.time()
        matrice=[]
        liste_espece=[]
        if phylum==None: #If Phylum is not specified we parse all existing phylums in the dataset
                for phylum in list(dictGeneral.keys()):
                        liste_espece+=list(dictGeneral[phylum].keys())
                        num=1
                        for i in dictGeneral[phylum]:
                                print(phylum+" Genome K-mer Computation :",num,"/",len(dictGeneral[phylum]))
                                pathTest=(dictGeneral[phylum][i])
                                dicSeq=KmerSignature(pathTest,kmer,normalized=True)
                                matrice.append(dicSeq)
                                num+=1
                matrice_Distance=np.zeros((len(matrice),len(matrice)))
                for i in range(len(matrice)):
                        for j in range(i,len(matrice)):
                                if i!=j:
                                        a=matrice[i].values[0]
                                        b=matrice[j].values[0]
                                        dst = distance.euclidean(a, b)
                                        matrice_Distance[i][j]=dst
                                        matrice_Distance[j][i]=dst

                matrice_distance_df=pd.DataFrame(data=matrice_Distance,columns=liste_espece,index=liste_espece)
                return matrice_distance_df
        else:
                liste_espece+=list(dictGeneral[phylum].keys())
                num=1
                for i in dictGeneral[phylum]:
                        print(phylum+" Genome K-mer Computation :",num,"/",len(dictGeneral[phylum]))
                        pathTest=(dictGeneral[phylum][i])
                        dicSeq=KmerSignature(pathTest,kmer,normalized=True)
                        matrice.append(dicSeq)
                        num+=1
                matrice_Distance=np.zeros((len(matrice),len(matrice)))
                for i in range(len(matrice)):
                        for j in range(i,len(matrice)):
                                if i!=j:
                                        a=matrice[i].values[0]
                                        b=matrice[j].values[0]
                                        dst = distance.euclidean(a, b)
                                        matrice_Distance[i][j]=dst
                                        matrice_Distance[j][i]=dst

                matrice_distance_df=pd.DataFrame(data=matrice_Distance,columns=liste_espece,index=liste_espece)
                return matrice_distance_df
        


"""Calcul Neighbour Joining"""

#Fonction récursive qui renvoie une séquence en format Newick
def Create_Tree(noeud):
    if noeud.fg=="" and noeud.fd=="":
        return(noeud.nom+":"+str(noeud.hauteur))
    else:
        arbre="("+"("+Create_Tree(noeud.fd)+","+Create_Tree(noeud.fg)+")"+":"+str(noeud.hauteur)+")"
        
    return arbre
    
#Mise à jour de la matrice et de la liste des noeuds pour chaque itération dans la fonction NJ
def Update_NJ_Matrix(matrice,listenoeud):
    liste_U=Sum_Distance_List(matrice)#Définition de l'ensemble des distances
    matrice_Q=Initialise_NJ_Matrix(matrice,liste_U)#Création de la matrice_Q de distance
    i,j=Compute_NJ_Minimal_Distance(matrice_Q)#Détermination des éléments les plus proches et les plus éloignés des autres éléments
    x,y=listenoeud[i].nom,listenoeud[j].nom
    noeudt=Create_NJ_Node(i,j,listenoeud,liste_U,matrice)
    matrice2=np.zeros((len(matrice)-1,len(matrice)-1))#Matrice qui permet la mise à jour
    
    t=0
    tt=0
    listenom=[]
    
    for it in listenoeud:
        listenom.append(it.nom)
    #Recopiage de la matrice de base dans la matrice à mettre à jour sans les éléments les plus proches
    for l in listenom:
        for c in listenom:
            matrice2[t][tt]=matrice.loc[l][c]
            tt=tt+1
        t=t+1
        tt=0
    #Ajout du nouveau noeud formé
    listenom.append(noeudt.nom)
    listenoeud.append(noeudt)
    matrice2_df=pd.DataFrame(matrice2,listenom,listenom)

    #Ajout des distances au nouveau noeud dans la matrice à mettre à jour
    for t in range(len(listenom)-1):
        matrice2_df.loc[noeudt.nom][listenom[t]]=(matrice.loc[x][listenom[t]]+matrice.loc[y][listenom[t]]-matrice.loc[x][y])/2
        matrice2_df.loc[listenom[t]][noeudt.nom]=(matrice.loc[x][listenom[t]]+matrice.loc[y][listenom[t]]-matrice.loc[x][y])/2
    return(matrice2_df)
    
#Création à partir d'une liste de noms des élements la liste des feuilles à la classe voulue 
def Create_NJ_Leaf(listenom):
    listefeuilles=[]
    for i in listenom:
        listefeuilles.append(NJ_Node(i,"","",0))
    return(listefeuilles)

#Création d'un nouveau noeud qui va remplacer les 2 valeurs les plus proches 
def Create_NJ_Node(i,j,listenoeud,liste_U,matrice):
    noeudt=NJ_Node(listenoeud[i].nom+listenoeud[j].nom,
          listenoeud[i],
          listenoeud[j],
          0)
    listenoeud[i].hauteur=(matrice.iloc[i][j]+liste_U[i]-liste_U[j])/2
    listenoeud[j].hauteur=(matrice.iloc[i][j]+liste_U[j]-liste_U[i])/2
    del listenoeud[i]
    del listenoeud[j]
    return noeudt
    
#Initialisation de la matrice Q qui pour chaque paire d'éléments renvoit une distance 
def Initialise_NJ_Matrix(matrice,liste_U):
    matriceQ=np.zeros((len(matrice),len(matrice)))
    for i in range(len(matrice)):
        for j in range(len(matrice)):
            if i!=j:
                matriceQ[i][j]=((matrice.iloc[i][j]-liste_U[i]-liste_U[j]))
    return(matriceQ) 
    
class NJ_Node:
    """Class qui défini un noeud tel que:
    -son nom
    -fils gauche
    -fils droit
    -sa hauteur
    -nombre d'élement gauche
    -nombre d'élement droit"""
    def __init__(self,nom,fg,fd,hauteur):
        self.nom=nom
        self.fg=fg
        self.fd=fd
        self.hauteur=hauteur

def Compute_NJ_Minimal_Distance(matrice):
    imin=1
    jmin=0
    minimal=matrice[imin][jmin]
    for i in range(1,len(matrice)):
        for j in range(0,i):
            if matrice[i][j]<minimal:
                minimal=matrice[i][j]
                imin=i
                jmin=j
    return(imin,jmin)

def Find_Column_Name(matrice):
    listenom=[]
    for i in range(1,len(matrice)+1):
        listenom.append(chr(i+64))
    return(listenom)
#Création de la liste qui stocke l'ensemble des distances de chaque élement
def Sum_Distance_List(matrice):
    liste_U=[]
    somme=0
    for i in range(len(matrice)):
        for j in range(len(matrice)):
            if i!=j:
                somme+=matrice.iloc[i][j]
        liste_U.append(somme/(len(matrice)-2))
        somme=0
    return(liste_U)

def NeighbourJoining(matrice_df:pd.DataFrame) -> str:
    """calculates Neighbor joining from distance matrix"""
    listenom=list(matrice_df.index)
    listenoeud=Create_NJ_Leaf(listenom)
    matrice_maj=Update_NJ_Matrix(matrice_df,listenoeud)
    #Si il ne reste qu'une valeur on ne peut pas la comparer à une autre donc on s'arrête quand il en reste 2
    while len(matrice_maj)>2:
        matrice_maj=Update_NJ_Matrix(matrice_maj,listenoeud)
    arbre=""
    for noeud in listenoeud:
        arbre+=Create_Tree(noeud)
    return arbre


"""Chi2 computation"""


def SequenceHomogeneity(path:str,kmer:int,fragmentSize:int)-> Tuple[list,list]:
        """Compute the Homogeneity of a sequence and spots any horizontal transfers """
        dicKmerSignature=KmerSignature(path,kmer,True)
        dicKmerSignature=dicKmerSignature/np.sum(dicKmerSignature.values)
        sequence=Read_Sequence(path)
        listepval = []
        listepos= []
        pos=0
        while pos<len(sequence):
                listepos.append(pos)
                print((pos/len(sequence))*100,"%")
                sequenceFragment=sequence[pos:pos+fragmentSize]
                seqCut = [sequenceFragment[i:i+kmer] for i in range(len(sequenceFragment)-(kmer-1)) ]
                dicSequenceHomogeneity = Count_Cuts(seqCut,False)
                dicComparaison = dicKmerSignature*np.sum(dicSequenceHomogeneity.values)
                resultat,pval = chisquare(np.array(dicSequenceHomogeneity.loc[0]),np.array(dicComparaison.loc[0]))
                listepval.append(pval)


                if pos+fragmentSize>len(sequence):
                        fragmentSize=len(sequence)-pos
                        pos+=fragmentSize
                else:
                        pos+=fragmentSize

        return listepval,listepos



if __name__=="__main__":
        #Test downloading sequences
        #DownloadSequences('./Bacteria.list','./Archea.list','./')
        
        #Test Parsing sequences
        #dictio=ParseSequences('./refseq/')
        
        #Test KmerSignature
        #print(KmerSignature('./testGenome.fna',3,True))
        
        #Test DistanceMatrix
        #testMatrix=DistanceMatrix(dictio,3)
        #print(testMatrix)
        
        #testMatrixPhylum=DistanceMatrix(dictio,3,'bacteria')
        #print(testMatrixPhylum)
        
        #Test NeighbourJoining
        #print(NeighbourJoining(testMatrix))
        
        #Test Chi2 Computation
        #print(SequenceHomogeneity('./testGenome.fna',3,1000))

