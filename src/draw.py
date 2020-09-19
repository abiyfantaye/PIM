import matplotlib.pyplot as plt    
import matplotlib.patches as patches
import matplotlib as mpl  
from matplotlib import gridspec
from scipy.interpolate import griddata
from scipy import interpolate
import numpy as np
import CWE as cwe

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
        oy = 0.75
        width = 0.75
        depth = 0.5

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        
        ax.plot([ox,ox], [0.1*oy, oy-depth/2.0], 'k-.', linewidth=1.25)
        border = patches.Rectangle((ox-width/2.0, oy-depth/2.0 ), width , depth, linewidth=linewidth,edgecolor='k',facecolor='grey')
        ax.add_patch(border)

        style="Simple,tail_width=0.75,head_width=6,head_length=12"
        kw = dict(arrowstyle=style, color="k")     
        
        arw1_r = oy-depth/2.0
        
        arw1_x1 = ox + arw1_r*np.sin(np.deg2rad(30))
        arw1_x2 = ox
        arw1_y1 = arw1_r - arw1_r*np.cos(np.deg2rad(30))
        arw1_y2 = oy-depth/2.0
        
        arw1 = patches.FancyArrowPatch(posA=(arw1_x1, arw1_y1), posB=(arw1_x2, arw1_y2), **kw)
        ax.text(0.1, 0.0, "Wind direction: $" + str(int(self.model.wind_direction)) + '^0$', fontsize=16)
        ax.add_patch(arw1)

        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        

    def _interpolate_cp(self, face, nx, ny, value_type, method='linear'):        

        cp = face.get_cp()
        x_tap, y_tap = face.get_coordinates() 
        
        
        x_grid = np.linspace(np.min(x_tap), np.max(x_tap), nx)
        y_grid = np.linspace(np.min(y_tap), np.max(y_tap), ny)
        
        
        extent = face.get_extent_coordinates()
        
        x_grid_face = np.linspace(extent[0,0], extent[1,0], nx)
        y_grid_face = np.linspace(extent[1,1], extent[2,1], nx)


        if value_type == 'mean':
            z = np.mean(cp, axis=1)
        if value_type == 'rms':
            z = np.std(cp, axis=1)
        if value_type == 'inst':
            z = cp[:,100]
            
        X, Y = np.meshgrid(x_grid, y_grid)
        Z = griddata((x_tap, y_tap), z, (X, Y), method=method)       
        f = interpolate.interp2d(x_grid, y_grid, Z, kind=method)
        
        
        X_face, Y_face = np.meshgrid(x_grid_face, y_grid_face)
        Z_face = f(x_grid_face, y_grid_face)
        return X_face, Y_face, Z_face
    
    def plot_wind_profile(self, value_type = 'mean'):
        #Mean velocity profile comparison with experement
        plt1, fig = cwe.setup_plot(plt, 22, 18, 18)
        
        ax = fig.add_subplot(1, 1, 1)    
        ax.tick_params(direction='in', size=8)    
        ax.set_xlim([0,1.25])
        ax.set_ylim([0,4])            
        H = self.model.building_height
        z = self.model.wind_profile[:,0]
        Uav = self.model.wind_profile[:,1]
        Iu = self.model.wind_profile[:,2]/100.0
        U_H = self.model.get_Uav(H)
        
        ax.set_ylabel('$z/H$')        
        ax.set_xlabel(r'$Uav/U_{H}, Iu$') 
        ax.plot(Uav/U_H, z/H,  'bo', markersize=8, markerfacecolor='none',markeredgewidth=1.5)
        ax.plot(Iu, z/H, 'rs', markersize=8, markerfacecolor='none',markeredgewidth=1.5)
            
        plt1.grid(linestyle='--')     
        
        fig.set_size_inches(15/2.54, 18/2.54)
        ax.legend(['Uav', '$Iu$'], loc=0) 
        plt1.tight_layout()
#        plt1.savefig('Plots/blwtl_vs_esdu.pdf')
#        plt1.savefig('Plots/blwtl_vs_esdu.png')
        plt1.show()
        
        
    def plot_cp(self, value_type = 'mean'):
        
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
        num_contours = 35
        
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
        
        markersize = 4
        
        gs = gridspec.GridSpec(1, n_plots, width_ratios=width_ratio, wspace=0.1, hspace=0.001) 
        gs0 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[4], height_ratios=[1.0, 3.0, 1.5])
        
        nx = 100
        ny = 100
        
        color_map = 'jet'
                
        if value_type == 'mean':
            data_range = np.linspace(-1.5, 1.0, num=num_contours)
            ticks=np.linspace(-1.5, 1.0, num=6)
        if value_type == 'rms':
            data_range = np.linspace(0.0, 0.5, num=num_contours)
            ticks=np.linspace(0.0, 0.6, num=5)
        if value_type == 'inst':
            data_range = np.linspace(-1.5, 1.5, num=num_contours)
            ticks=np.linspace(-1.5, 1.0, num=6)

        interp_method = 'linear'

        #plot the four side faces
        for i in range(n_plots-1):   
            face = self.model.faces[i]         
            ax = fig.add_subplot(gs[i])
            x,y = face.get_coordinates()
            ax.plot(scale*x, scale*y, 'k+', markersize=markersize)
            
            X,Y,Z = self._interpolate_cp(face, nx, ny, value_type, method=interp_method)            
            
#            cs = ax.contour(scale*X, scale*Y, Z, data_range[0:num_contours:4], colors='k',alpha=0.75, linewidths=0.75, linestyles='solid')            
            cs = ax.contour(scale*X, scale*Y, Z, data_range, colors='k',alpha=0.75, linewidths=0.75, linestyles='solid')            
#            ax.clabel(cs, fmt='%2.1f', colors='k', fontsize=14)

            cs = ax.contourf(scale*X, scale*Y, Z, data_range, cmap=color_map)

            ax.set_title(face.plot_name)
            extent = face.get_extent_coordinates()
            ax.set_xlim(scale*extent[0,0], scale*extent[1,0])
            ax.set_ylim(scale*extent[1,1], scale*extent[2,1])  
            if face.name != 'North':
                ax.set_yticks([])            
                
        #plot the roof 
        if self.plot_top == True:
            face = self.model.faces[-1]             
            ax = fig.add_subplot(gs0[0])
            x,y = face.get_coordinates()
            ax.plot(scale*x, scale*y, 'k+', markersize=markersize)
            
            X,Y,Z = self._interpolate_cp(face, nx, ny, value_type, method=interp_method)
            cs = ax.contour(scale*X, scale*Y, Z, data_range, colors='k',alpha=0.75, linewidths=0.75, linestyles='solid')            

            cs = ax.contourf(scale*X, scale*Y, Z, data_range, cmap=color_map)

            ax.set_title(face.plot_name)
            extent = face.get_extent_coordinates()
            ax.set_xlim(scale*extent[0,0], scale*extent[1,0])
            ax.set_ylim(scale*extent[1,1], scale*extent[2,1])  
            ax.yaxis.set_ticks_position('right')            
            
        #plot the legend 
        ax = fig.add_subplot(gs0[1])
        ax.set_visible(False)
        cbar = fig.colorbar(cs, ax=ax, fraction=1.0, aspect=12.5, ticks=ticks)
        cbar.ax.set_ylabel('$C_p$')   
        
        #plot the key 
        ax = fig.add_subplot(gs0[2])
        self.plot_key(ax)
        
        
        fig.set_size_inches(35/2.54, 50/2.54)
        plt.tight_layout()
        plt.show()