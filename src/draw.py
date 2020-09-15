import matplotlib.pyplot as plt    
import matplotlib.patches as patches
import matplotlib as mpl  
from matplotlib import gridspec

import numpy as np

class Plotter:

    def __init__(self, model, plot_top=True):
        self.model = model 
        self.plot_top = plot_top
        
    def _plot_face_frame(self, ax, face):
        ax.axis('off')
        border = patches.Rectangle(-face.width/2.0,0.0 ,face.width , face.height)
        ax.add_patch(border)
        ax.set_title(face.name)
        ax.set_ylim(0.0, face.height)   
    
    def plot(self):
        

        n_plots = len(self.model.faces)
        if self.plot_top == False:
            n_plots -= 1
            
        fig = plt.figure(facecolor='white')
        markersize = 5
        linewidth = 1.25   
        
        #ratio to adjust proportion of each subplot
        width_ratio = np.zeros(n_plots)

        for i in range(n_plots):        
            width_ratio[i] = self.model.faces[i].width/self.model.faces[0].width
        
        gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratio, wspace=0.1) 
        gs0 = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[4])

        for i in range(n_plots):  
            if i == n_plots-1:
                ax = fig.add_subplot(gs0[0])
            else :
                ax = fig.add_subplot(gs[i])
            face = self.model.faces[i]            
            ax.set_title(face.name)
            ax.set_xlim(0.00,face.width) 
            ax.set_ylim(0.00,face.height) 
            ax.set_xticks([])
            ax.set_yticks([])

        plt.tight_layout()
        plt.show()
#        
#        
#        for i in range(len(self.north)):
#            ax.plot(self.taps[self.north[i]].y, self.taps[self.north[i]].z, 'k+', markersize=markersize)
##            ax.text(self.taps[self.north[i]].y, self.taps[self.north[i]].z, self.taps[self.north[i]].tag, fontsize=7,rotation=45)
#            
#        ax = fig.add_subplot(1, n_plots, 2)            
#        ax.axis('off')
#        border = patches.Rectangle((-self.depth/2.0,0.0 ), self.depth , self.height, linewidth=linewidth,edgecolor='k',facecolor='none')
#        ax.add_patch(border)
#        ax.set_title('West')
#        ax.set_ylim(0.00,self.height) 
#        for i in range(len(self.west)):
#            ax.plot(self.taps[self.west[i]].x, self.taps[self.west[i]].z, 'k+', markersize=markersize)
##            ax.text(self.taps[self.west[i]].x, self.taps[self.west[i]].z, self.taps[self.west[i]].tag, fontsize=7,rotation=45)
##
##            
#        ax = fig.add_subplot(1, n_plots, 3)            
#        ax.axis('off')
#        border = patches.Rectangle((-self.width/2.0,0.0 ), self.width , self.height, linewidth=linewidth,edgecolor='k',facecolor='none')
#        ax.add_patch(border)
#        ax.set_title('South')            
#        ax.set_ylim(0.00,self.height)   
#        for i in range(len(self.south)):
#            ax.plot(self.taps[self.south[i]].y, self.taps[self.south[i]].z, 'k+', markersize=markersize)
##            ax.text(self.taps[self.south[i]].y, self.taps[self.south[i]].z, self.taps[self.south[i]].tag, fontsize=7,rotation=45)            
##
#        ax = fig.add_subplot(1, n_plots, 4)           
#        ax.axis('off')
#        border = patches.Rectangle((-self.depth/2.0,0.0 ), self.depth , self.height, linewidth=linewidth,edgecolor='k',facecolor='none')
#        ax.add_patch(border)
#        for i in range(len(self.east)):
#            ax.plot(self.taps[self.east[i]].x, self.taps[self.east[i]].z, 'k+', markersize=markersize)
##            ax.text(self.taps[self.east[i]].x, self.taps[self.east[i]].z, self.taps[self.east[i]].tag, fontsize=7,rotation=45)
##
##
##        ax = fig.add_subplot(1, n_plots, 5)            
##        ax.axis('off')
##        border = patches.Rectangle((-self.depth/2.0, -self.width/2.0), self.depth , self.width, linewidth=linewidth, edgecolor='r', facecolor='none')
##        ax.add_patch(border)        
##        ax.set_title('East')            
##        for i in range(len(self.top)):
##            ax.plot(self.taps[self.top[i]].x, self.taps[self.top[i]].y, 'k+', markersize=markersize)
##            ax.text(self.taps[self.top[i]].x, self.taps[self.top[i]].y, self.taps[self.top[i]].tag, fontsize=7,rotation=90)
#        
#        fig.set_size_inches(30/2.54, 75/2.54)    
#        plt.tight_layout()
#        plt.show()
#
# 