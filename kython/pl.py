import matplotlib.pyplot as plt
from matplotlib.patches import *
import numpy as np
import pandas as pd
from ete3 import Tree, TreeStyle, AttrFace, faces, NodeStyle
import os




def KmerSignature(signature_df:pd.DataFrame,organism:str)->None:
  """
  signature_df: DataFrame or list of Dataframe produced by KmerSignature
  Plotting function of Kmer signature, will do multiple plot if multiple signatures in list
  """

  if type(signature_df)==list:
    for i in range(len(signature_df)):
      # Compute pie slices
      N = len(signature_df[i].columns)
      theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
      radii = signature_df[i].loc[0].to_numpy()
      width = 2* np.pi / N
      colors = plt.cm.viridis(100)
      ax_tot = len(signature_df)
      ax = plt.subplot(100+ax_tot*10+i+1, projection='polar')
      ax.bar(theta, radii, width=width, bottom=0.0, color=colors, alpha=0.8)
      ax.set_yticklabels([])
      ax.set_xticklabels([])
  else:
    N = len(signature_df.columns)
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    radii = signature_df.loc[0].to_numpy()
    width = 2* np.pi / N
    colors = plt.cm.viridis(100)
    ax = plt.subplot(111, projection='polar')
    ax.bar(theta, radii, width=width, bottom=0.0, color=colors, alpha=0.8)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
  ax.set_title(organism+" Kmer Signature")
  print()
  ax.text(0.030,0.045,"Kmer size: "+str(len(signature_df.columns.values[0]))+" pb" )


  plt.show()
  
  

def SequenceHomogeneity(pValueList:list,posList:list,pValueAccepted:float,organism:str,kmerSize:int=None) -> None:
        


        #Starting points for shapes
        startingX=0
        startingY=5
        
        #Shapes configuration 
        blueRectangleHeight=2
        blueRectangleLength=100
        blueRectangleColor=plt.cm.viridis(100)
        
        redRectanglesHeight=2
        redRectanglesColor='#f03434'
        
        ticksColor="#17181a"
        ticksHeight=2
        ticksThickness=0.5
        
        lengthGenome=posList[-1]/100
        
        #Creating shapes
        
        #Big Rectangle
        blueRectangle=Rectangle((startingX,startingY), width=blueRectangleLength,height=blueRectangleHeight,color=blueRectangleColor)
        
        #Tiny Rectangles
        listRedRectangles=[]
        for i in range(len(pValueList)):
                if pValueList[i]<pValueAccepted:
                        listRedRectangles.append(Rectangle((posList[i]/lengthGenome,startingY), width=(posList[i+1]-posList[i])/lengthGenome,height=redRectanglesHeight,color=redRectanglesColor))
                        
        
        #Ticks
        listTicks=[]
        
        for positionScaled in np.linspace(0,100,10):
                line=plt.Line2D((positionScaled,positionScaled),(startingY,startingY+ticksHeight),linewidth=ticksThickness,color=ticksColor)
                listTicks.append(line)
                
                

        
        
        
        
        
        #Creating the figure
        fig=plt.figure(figsize=(20,4))#figsize=(blueRectangleLength*2,blueRectangleHeight*2))
        ax=fig.add_subplot(111)
        
        #Annotations
        listAnnotationsPos=np.linspace(0,100,4)
        listAnnotations=[]
        for annotation in listAnnotationsPos:
                ax.text(annotation,3.5,str(int(annotation * lengthGenome)),color=ticksColor,ha='center')
                
        ax.text(50,2,"Genome Length (pb)",color=ticksColor,ha='center')
        
        ax.text(75,0,"Percentage of Horizontal transfers: "+str(format(len(listRedRectangles)/(len(posList)-1),'.2f')),color='black')
        if kmerSize:
                ax.text(25,0,"Kmer size: "+str(kmerSize),color='black')
                
        #Legend
        legend_elements=[Patch(color=blueRectangleColor,label="Corresponding Kmer Signature"),
                         Patch(color=redRectanglesColor,label= "Horizontal transfert")]
        ax.legend(handles=legend_elements,loc='upper right',bbox_to_anchor=(1,4))
        
        #Title
        ax.set_title(organism,color=ticksColor)
        

        #Adding shapes to the ax
        ax.add_patch(blueRectangle)
        
        for tick in listTicks:
                ax.add_line(tick)
                
        for redRectangle in listRedRectangles:
                ax.add_patch(redRectangle)
        
        #Adjusting ax
        ax.axis("scaled")
        ax.get_yaxis().set_visible(False)
        ax.axis("off")
        
        
        plt.show()


def GraphTree(treeNewick:str, pathDictionnary:dict, outputPath:str) -> None:
  """
  treeNewick: Newick tree format (from NJ function)
  pathDictionnary: To get the class corresponding to name for coloring
  outputPath: name to give to the .png tree figure produced

  Coloring: Red = archaea;  Bleu = Bacteria
  """

  os.environ['QT_QPA_PLATFORM']='offscreen'

  keys_archaea = []
  for i in list(pathDictionnary['archaea'].keys()):
    keys_archaea.append(i.replace(' ','_'))

  treeNewick = treeNewick.replace(')((',',(')+'OROOT;'
  t = Tree(treeNewick, format=1)

  ts = TreeStyle()
  ts.mode = "c"

  style2 = NodeStyle()
  style2["fgcolor"] = "#2e3131"
  style2["vt_line_color"] = "#2e3131"
  style2["hz_line_color"] = "#2e3131"
  style2["vt_line_width"] = 20
  style2["hz_line_width"] = 20
  style2["vt_line_type"] = 0 
  style2["hz_line_type"] = 0

  def layout(node):
    """
    Hidden function to build the layout of the Tree. Make the coloring of the 
    leaves according to the class bacteria/archaea
    """
    if node.is_leaf():
        if node.name in keys_archaea:
          N = AttrFace("name", ftype='Arial', fsize=30, fgcolor='black')
          N.background.color='Crimson'

        else:
          N = AttrFace("name", ftype='Arial', fsize=30, fgcolor='black', penwidth=2)
          N.background.color='CornflowerBlue'

        faces.add_face_to_node(N, node, 0, aligned=True, position="aligned")
        node.img_style = style2


  ts.layout_fn = layout
  ts.show_leaf_name = False

  t.render(outputPath+'.png', w=1000,tree_style=ts)



if __name__=='__main__':

        import pp

        #Test KmerSignature
        #signature=pp.KmerSignature('./testGenome.fna',3,True)
        #KmerSignature(signature,'Poulet')

        #Test Sequence Homogeneity
        #listpVal,listPos=pp.SequenceHomogeneity('./testGenome.fna',5,30000)
        #SequenceHomogeneity(listpVal,listPos,0.05,'Poulet',3)
        
        #Test GraphTree
        dictio=pp.ParseSequences('./refseq/')
        print(dictio)
        
        testMatrix=pp.DistanceMatrix(dictio,3)
        tree=pp.NeighbourJoining(testMatrix)
        GraphTree(tree, dictio,'lolilol.png')
        
        pass
