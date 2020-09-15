import matplotlib.pyplot as plt    
import matplotlib.patches as patches
import matplotlib as mpl  
from matplotlib import gridspec

import numpy as np

class Plotter:

    def __init__(self, model, plot_top=True):
        self.model = model 
        self.plot_top = plot_top
        
    def _get_face_coordinates(self, face):
        n_taps  = len(face.taps)
        x = np.zeros(n_taps)
        y = np.zeros(n_taps)
        
        for i in range(n_taps):
            if face.name == 'North':
                x[i] = -face.taps[i].coord.y 
                y[i] = face.taps[i].coord.z 
            if face.name == 'West':
                x[i] = face.taps[i].coord.x  
                y[i] = face.taps[i].coord.z 
            if face.name == 'South':
                x[i] = face.taps[i].coord.y 
                y[i] = face.taps[i].coord.z 
            if face.name == 'East':
                x[i] = -face.taps[i].coord.x 
                y[i] = face.taps[i].coord.z 
            if face.name == 'Top':
                x[i] = face.taps[i].coord.y 
                y[i] = -face.taps[i].coord.x
                        
        return x, y

    
    def plot(self):
        
        n_plots = len(self.model.faces)
        scale = self.model.scale
        
        #ratio to adjust proportion of each subplot
        width_ratio = np.zeros(n_plots)

        for i in range(n_plots):        
            width_ratio[i] = self.model.faces[i].width/self.model.faces[0].width
                
        fig = plt.figure(facecolor='white')
        markersize = 5
        linewidth = 1.25  
        
        gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratio, wspace=0.1, hspace=0.001) 
        gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[4], height_ratios=[1.0, 6.0])

        #plot the four side faces
        for i in range(n_plots-1):   
            face = self.model.faces[i]         
            ax = fig.add_subplot(gs[i])
            x,y = self._get_face_coordinates(face)
            ax.plot(scale*x, scale*y, '+', markersize=markersize)
            
#            for j in range(len(face.taps)):
#                ax.text(scale*x[j], scale*y[j], face.taps[j].name, fontsize=8,rotation=45)

            ax.set_title(face.name + ' Wall')
            
            ax.set_xlim(-scale*face.width/2.0, scale*face.width/2.0)
            ax.set_ylim(0.00,scale*face.height)  
            
            if face.name != 'North':
                ax.set_yticks([])

        
        #plot the roof 
        if self.plot_top == True:
            face = self.model.faces[-1]             
            ax = fig.add_subplot(gs0[0])
            x,y = self._get_face_coordinates(face)
            ax.plot(scale*x, scale*y, '+', markersize=markersize)
                
    #            for j in range(len(face.taps)):
    #                ax.text(scale*x[j], scale*y[j], face.taps[j].name, fontsize=8,rotation=45)
    
            ax.set_title('Roof')
            ax.set_xlim(-scale*face.width/2.0, scale*face.width/2.0)
            ax.set_ylim(-scale*face.height/2.0, scale*face.height/2.0)

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