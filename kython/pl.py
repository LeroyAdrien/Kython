import matplotlib.pyplot as plt
from matplotlib.patches import *
import pp

def SequenceHomogeneity(pValueList:list,posList:list,pValueAccepted:float) -> None:
        #Creating figure
        fig, ax = plt.subplots()
        
        #Shapes configuration 
        blueRectangleHeight=0.3
        blueRectangleLength=100
        blueRectangleColor='blue'
        
        redRectanglesHeight=0.3
        redRectanglesColor='red'
        
        #Creating shapes
        blueRectangle=Rectangle((0,0), width=blueRectangleLength,height=blueRectangleHeight,color=blueRectangleColor)
        
        #Adding shapes to the ax
        ax.add_patch(blueRectangle)

        fig.show()



if __name__=='__main__':
        
        listpVal,listPos=pp.SequenceHomogeneity('./testGenome.fna',3,1000)
        SequenceHomogeneity(listpVal,listPos,0.05)
