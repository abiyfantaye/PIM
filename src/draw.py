import matplotlib.pyplot as plt    
import matplotlib.patches as patches
import matplotlib as mpl  
from matplotlib import gridspec
from scipy.interpolate import griddata
from scipy import interpolate
import numpy as np

class Plotter:

    def __init__(self, model, plot_top=True):
        self.model = model 
        self.plot_top = plot_top

    def plot_key(self, ax):
        """
        Plots the key showing the building and wind direction.  
        """
        linewidth = 1.5
        ox = 0.5
        oy = 0.25
        width = 0.8
        depth = 0.15

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        ax.plot([ox,ox], [0.15*oy, oy-depth/2.0], 'k-.', linewidth=1.25)
        border = patches.Rectangle((ox-width/2.0, oy-depth/2.0 ), width , depth, linewidth=linewidth,edgecolor='k',facecolor='grey')
        ax.add_patch(border)

        style="Simple,tail_width=0.75,head_width=6,head_length=12"
        kw = dict(arrowstyle=style, color="k")     
        
        arw1_r = oy-depth/2.0
        
        arw1_x1 = ox + 1.5*arw1_r*np.sin(np.deg2rad(45))
        arw1_x2 = ox
        arw1_y1 = arw1_r - 1.1*arw1_r*np.cos(np.deg2rad(45))
        arw1_y2 = oy-depth/2.0
        
        arw1 = patches.FancyArrowPatch(posA=(arw1_x1, arw1_y1), posB=(arw1_x2, arw1_y2), **kw)
        ax.text(0.1, 0.0, "Wind direction: $" + str(int(self.model.wind_direction)) + '^0$', fontsize=16)
        ax.add_patch(arw1)

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)

    def _interpolate_cp(self, face, nx, ny):        

        cp = face.get_cp()
        x_tap, y_tap = face.get_coordinates() 
        
        x_grid = np.linspace(np.min(x_tap), np.max(x_tap), nx)
        y_grid = np.linspace(np.min(y_tap), np.max(y_tap), ny)
        
        x_grid_face = np.linspace(-face.width/2.0, face.width/2.0, nx)
        
        if face.name == "Top":
            y_grid_face = np.linspace(-face.height/2.0, face.height/2.0, nx)
        else:
            y_grid_face = np.linspace(0.0, face.height, ny)

        z = np.mean(cp, axis=1)    
        X, Y = np.meshgrid(x_grid, y_grid)
        Z = griddata((x_tap, y_tap), z, (X, Y), method='linear')       
        f = interpolate.interp2d(x_grid, y_grid, Z, kind='linear')
        
        
        X_face, Y_face = np.meshgrid(x_grid_face, y_grid_face)
        Z_face = f(x_grid_face, y_grid_face)
        return X_face, Y_face, Z_face
        
    def plot(self):
        
        n_plots = len(self.model.faces)
        scale = self.model.scale
        
        #ratio to adjust proportion of each subplot
        width_ratio = np.zeros(n_plots)

        for i in range(n_plots):        
            width_ratio[i] = self.model.faces[i].width/self.model.faces[0].width
                
        fig = plt.figure(facecolor='white')
        
        font_size=18
        legend_font_size=20
        axis_font_size=14
        
        font = {'family' : 'Times New Roman','weight' : 'normal', 'size'   : font_size}
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        plt.rc('font', **font)    
        plt.rc('axes', labelsize=font_size)    # fontsize of the x and y labels
        plt.rc('axes', titlesize=font_size)  # fontsize of the axes title
        plt.rc('xtick', labelsize=axis_font_size)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=axis_font_size)    # fontsize of the tick labels
        plt.rc('axes', linewidth=1.25)    
        plt.rc('legend', fontsize=legend_font_size)
        plt.rc('text', usetex=True)
        
        markersize = 5
        linewidth = 1.25  
        
        gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratio, wspace=0.1, hspace=0.001) 
        gs0 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[4], height_ratios=[1.0, 5.0])
        
        nx = 100
        ny = 100
        
        color_map = 'jet'

        #plot the four side faces
        for i in range(n_plots-1):   
            face = self.model.faces[i]         
            ax = fig.add_subplot(gs[i])
            x,y = face.get_coordinates()
            ax.plot(scale*x, scale*y, '+', markersize=markersize)
            
            X,Y,Z = self._interpolate_cp(face, nx, ny)
            ax.contourf(scale*X, scale*Y, Z, cmap=color_map)

            ax.set_title(face.name + ' Wall')
            ax.set_xlim(-scale*face.width/2.0, scale*face.width/2.0)
            ax.set_ylim(0.00,scale*face.height)  
            if face.name != 'North':
                ax.set_yticks([])
            
            for j in range(len(face.taps)):
                ax.text(scale*x[j], scale*y[j], face.taps[j].name, fontsize=12,rotation=45)
                
        #plot the roof 
        if self.plot_top == True:
            face = self.model.faces[-1]             
            ax = fig.add_subplot(gs0[0])
            x,y = face.get_coordinates()
            ax.plot(scale*x, scale*y, '+', markersize=markersize)
            
            X,Y,Z = self._interpolate_cp(face, nx, ny)
            ax.contourf(scale*X, scale*Y, Z, cmap=color_map)

            ax.set_title('Roof')
            ax.set_xlim(-scale*face.width/2.0, scale*face.width/2.0)
            ax.set_ylim(-scale*face.height/2.0, scale*face.height/2.0)
            ax.yaxis.set_ticks_position('right')            
            
            for j in range(len(face.taps)):
                ax.text(scale*x[j], scale*y[j], face.taps[j].name, fontsize=12,rotation=45)
            
        #plot the key 
            
        ax = fig.add_subplot(gs0[1])
        self.plot_key(ax)

            
        fig.set_size_inches(35/2.54, 50/2.54)

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