import numpy as np
from scipy import stats
import pandas as pd
from scipy import signal
from tap import Tap
from point import Point 
from face import Face 
import CWE as cwe


   
class PIM:
    """A class representing Pressure Integration Model(PIM)"""

    def __init__(self, 
                 cp_file_name, 
                 tap_file_name,
                 wind_direction,
                 building_width, 
                 building_depth,                  
                 building_height, 
                 z0, u_ref, 
                 z_ref, 
                 gradient_height, 
                 gradient_wind_speed, 
                 scale, 
                 broken_taps):
        """Constructs PIM.
        Args:
          cp_file_name: name of the file that contains the Cp data, the column
              representing the taps and rows the time series.
          tap_file_name: name of the file that contains the Cp data, the column.
          .
          .
          .
          .
          """
        self.cp_file_name = cp_file_name
        self.tap_file_name = tap_file_name
        self.wind_direction = wind_direction
        self.building_height = building_height
        self.building_width = building_width
        self.building_depth = building_depth
        self.z0 = z0
        self.u_ref = u_ref
        self.z_ref = z_ref            
        self.gradient_height = gradient_height
        self.gradient_wind_speed = gradient_wind_speed
        self.scale = scale
        self.broken_taps = broken_taps
        self.faces = []
        self.taps  = []    
        self.__create_taps()
        self.__create_faces()
        
#        self.__read_cp_data()
#        self.__correct_cp_to_building_height()

    def __read_cp_data(self):
        cp_data = np.loadtxt(self.Cp_file_path)

        self.cp_data = np.zeros((self.tap_count, np.shape(cp_data)[1]))
        
        for i in range(self.tap_count):
            self.cp_data[i,:] = cp_data[i,:]
        
        #Take only the first 120second 
#        self.cp_data = self.cp_data[:,0:48000]

    def __create_taps(self):
        """
        Creates taps reading information from a text file. The tap coordinate should be formated as: 
        
        TapID   Face    X-coord     Y-Coord     Z-coord
        -----   -----   ------      -------     -------
        """
        tap_file = open(self.tap_file_name, "r")
        lines  = tap_file.readlines()
        tap_index = 0
        for line in lines:
            atribute = line.split('\t')
            coord  = Point(float(atribute[2]), float(atribute[3]),  float(atribute[4]))
            tap = Tap(tap_index, atribute[0], coord, atribute[1])
            self.taps.append(tap)
            tap_index += 1   
        
        self.tap_count = len(self.taps)

#        self.__assign_tap_faces()
#        self.__create_neighbouring_taps()
        
    def __create_faces(self):    
        """
        Creates each faces as 'North', 'South', 'East', 'West', 'Top'.
        
        """
        #Create each face
        
        self.faces.append(Face('North', self.building_width,self.building_height, Point(-1.0, 0.0, 0.0)))        
        self.faces.append(Face('West', self.building_depth,self.building_height, Point(0.0, -1.0, 0.0)))
        self.faces.append(Face('South', self.building_width,self.building_height, Point(1.0, 0.0, 0.0)))
        self.faces.append(Face('East', self.building_depth,self.building_height, Point(0.0, 1.0, 0.0)))
        self.faces.append(Face('Top', self.building_width, self.building_depth, Point(0.0, 0.0, 1.0)))
        
        self.face_count  = len(self.faces)

        #Assign each tap to corresponding faces
        for i in range(self.tap_count):
            for j in range(self.face_count):
                if self.taps[i].face == self.faces[j].name:
                    self.faces[j].taps.append(self.taps[i])

                
#    def __create_neighbouring_taps(self):    
#        """
#        Creates the neighbouring taps depending on the face where the tap is located.
#                
#        """
#        dist = 1.0e20
#        top = ''
#        bottom = ''
#        left = ''
#        right = ''
#        
#        for tapi in self.taps:
#            for tapj in self.taps:
#                temp_dist =  calculate_distances(tapi.x, tapi.y, tapi.z, tapj.x, tapj.y, tapj.z)
#                if temp_dist < dist:
#                    if tapi.face == tapj.face:
#                        if tapi.face == 'North':
#                            if tapi.y == tapj.y and  tapi.z < tapj.z:
#                                top = tapj.tag
#                            if tapi.y == tapj.y and  tapi.z > tapj.z:
#                                bottom = tapj.tag           
#                            if tapi.z == tapj.z and  tapi.y > tapj.y:
#                                right = tapj.tag          
#                            if tapi.z == tapj.z and  tapi.y < tapj.y:
#                                left = tapj.tag          
#                        if tapi.face == 'South':
#                            if tapi.y == tapj.y and  tapi.z < tapj.z:
#                                top = tapj.tag
#                            if tapi.y == tapj.y and  tapi.z > tapj.z:
#                                bottom = tapj.tag           
#                            if tapi.z == tapj.z and  tapi.y > tapj.y:
#                                right = tapj.tag          
#                            if tapi.z == tapj.z and  tapi.y < tapj.y:
#                                left = tapj.tag  
#                        if tapi.face == 'West':
#                            if tapi.x == tapj.x and  tapi.z < tapj.z:
#                                top = tapj.tag
#                            if tapi.x == tapj.x and  tapi.z > tapj.z:
#                                bottom = tapj.tag           
#                            if tapi.z == tapj.z and  tapi.x > tapj.x:
#                                right = tapj.tag          
#                            if tapi.z == tapj.z and  tapi.x < tapj.x:
#                                left = tapj.tag  
#                        if tapi.face == 'East':
#                            if tapi.x == tapj.x and  tapi.z < tapj.z:
#                                top = tapj.tag
#                            if tapi.x == tapj.x and  tapi.z > tapj.z:
#                                bottom = tapj.tag           
#                            if tapi.z == tapj.z and  tapi.x > tapj.x:
#                                right = tapj.tag          
#                            if tapi.z == tapj.z and  tapi.x < tapj.x:
#                                left = tapj.tag  
#                        if tapi.face == 'Top':
#                            if tapi.x == tapj.x and  tapi.y < tapj.y:
#                                top = tapj.tag
#                            if tapi.x == tapj.x and  tapi.y > tapj.y:
#                                bottom = tapj.tag           
#                            if tapi.y == tapj.y and  tapi.x > tapj.x:
#                                right = tapj.tag          
#                            if tapi.y == tapj.y and  tapi.x < tapj.x:
#                                left = tapj.tag  
#            tapi.top = top
#            tapi.bottom = bottom
#            tapi.left = left
#            tapi.right = right
#            tapi.neighbours = [top, bottom, left, right]
#            dist = 1.0e20
#            top = ''
#            bottom = ''
#            left = ''
#            right = ''            
#            
#            
#    def __get_face_cp(self, face):
#        
#        if face == 'North':
#            x_coord = np.zeros(len(self.north))
#            y_coord = np.zeros(len(self.north))
#            cp = np.zeros((len(self.north), np.shape(self.cp_data)[1]))
#            for i in range(len(self.north)):
#                x_coord[i] = self.taps[self.north[i]].y
#                y_coord[i] = self.taps[self.north[i]].z
#                cp[i,:] = self.cp_data[self.north[i],:]
#        if face == 'West':
#            x_coord = np.zeros(len(self.west))
#            y_coord = np.zeros(len(self.west))
#            cp = np.zeros((len(self.west), np.shape(self.cp_data)[1]))
#            for i in range(len(self.west)):
#                x_coord[i] = self.taps[self.west[i]].x
#                y_coord[i] = self.taps[self.west[i]].z
#                cp[i,:] = self.cp_data[self.west[i],:]
#        if face == 'South':
#            x_coord = np.zeros(len(self.south))
#            y_coord = np.zeros(len(self.south))
#            cp = np.zeros((len(self.south), np.shape(self.cp_data)[1]))
#            for i in range(len(self.south)):
#                x_coord[i] = self.taps[self.south[i]].y
#                y_coord[i] = self.taps[self.south[i]].z
#                cp[i,:] = self.cp_data[self.south[i],:]
#        if face == 'East':
#            x_coord = np.zeros(len(self.east))
#            y_coord = np.zeros(len(self.east))
#            cp = np.zeros((len(self.east), np.shape(self.cp_data)[1]))
#            for i in range(len(self.east)):
#                x_coord[i] = self.taps[self.east[i]].x
#                y_coord[i] = self.taps[self.east[i]].z
#                cp[i,:] = self.cp_data[self.east[i],:]                     
#        
#        return x_coord, y_coord, cp
#    
#    def find_tap_by_tag(self, tag):
#        found_tap = -1 
#        
#        for tap in self.taps:
#            if tap.tag == tag:
#                found_tap = tap.idx
#                break
#        
#        return found_tap
#        
#    def __correct_cp_to_building_height(self):
#        
#        #Correct the velocity using the velocity ratio of 
#        #ESDU profile at two points, the gradient 
#        uref = 20.0
#        zref = 10.0
#
#        esdu = ESDU(uref, zref, self.z0)
#        z_gradient = 1.4732 # Gradient wind tunnle height
#        u_grdient = esdu.get_V_z(self.scale*z_gradient)
#        u_h = esdu.get_V_z(self.scale*self.height)
#        
#        corr = (u_grdient/u_h)**2.0
#        
#        print(corr)
#        
#        self.cp_data = self.cp_data*corr
#        
#    def __interpolate_tap_data(self, ngridx, ngridy, face):        
#        from scipy.interpolate import griddata
#        from scipy import interpolate
#        x_coord, y_coord, cp = self.__get_face_cp(face)
#        x_grid = np.linspace(np.min(x_coord), np.max(x_coord), ngridx)
#        y_grid = np.linspace(np.min(y_coord), np.max(y_coord), ngridy)
#        
#        x_grid_new = np.linspace(-self.width/2.0, self.width/2.0, ngridx)
#        y_grid_new = np.linspace(0.0, self.height, ngridy)
#        z = np.mean(cp, axis=1)    
#        X, Y = np.meshgrid(x_grid, y_grid)
#        Z = griddata((x_coord, y_coord), z, (X, Y), method='linear')       
#        f = interpolate.interp2d(x_grid, y_grid, Z, kind='linear')
#        X_new, Y_new = np.meshgrid(x_grid_new, y_grid_new)
#        Z_new = f(x_grid_new, y_grid_new)
#        return X_new, Y_new, Z_new
#    
##    def __correct_broken_taps(self):
#        
#        #If some of the taps are broken, this fucntion corrects the 
#        #pressure data by taking average between the nearest taps
#        
##        for i in range(len(self.broken_taps)):
##            tap = self.__find_tap_by_tag(self.broken_taps[i])
#
#        
##    def __find_neighbouring_taps(self, tap_tag):
#        
##        tap = self.__find_tap_by_tag(tap_tag)
#        
##        if tap.face == 'North': 
#        
#    def plot_faces_with_taps(self):
#        import matplotlib.pyplot as plt    
#        import matplotlib.patches as patches
#        import matplotlib as mpl    
#           
#        fig = plt.figure(facecolor='white')
#        markersize = 5
#        linewidth = 1.25
#        n_plots = 5
#        
#        ax = fig.add_subplot(1, n_plots, 1)
#
#        ax.axis('off')
#        border = patches.Rectangle((-self.width/2.0,0.0 ),self.width , self.height, linewidth=linewidth,edgecolor='k',facecolor='none')
#        ax.add_patch(border)
#        ax.set_title('North')
#        ax.set_ylim(0.00, self.height)   
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
#    def plot_taps_and_walls(self):
#        import matplotlib.pyplot as plt    
#        import matplotlib.patches as patches
#        import matplotlib as mpl    
#
##        tap1 = self.__find_tap_by_tag('609')
##        tap2 = self.__find_tap_by_tag('608')
##        
##        temp = self.cp_data[tap1.idx,:]
##        
##        self.cp_data[tap1.idx,:] = self.cp_data[tap2.idx,:]
##        self.cp_data[tap2.idx,:] = temp 
#
##        plt.plot(self.cp_data[tap1.idx,:])
##        plt.plot(self.cp_data[tap2.idx,:])
#           
#        fig = plt.figure(facecolor='white')
#        markersize = 3
#        linewidth = 0.25
#        n_plots = 4
#        ngridx = 70
#        ngridy = 100
#        ncontour = 20
#        vmin = 0.0
#        vmax = 1.0
#        
#        ax = fig.add_subplot(1, n_plots, 1)
#        xx, yy, zz = self.__interpolate_tap_data(ngridx, ngridy, 'North')
#        
##        ax.pcolorfast(xx, yy, zz, cmap='RdBu_r', vmin=vmin, vmax=vmax)
#        
#        ax.contour(xx, yy, zz, ncontour, linewidths=0.5, colors='k', vmin=vmin, vmax=vmax)
#        ctrplt = ax.contourf(xx, yy, zz, ncontour, cmap="RdBu_r", vmin=vmin, vmax=vmax)        
#        plt.colorbar(ctrplt, ax=ax)
#
#
#        for i in range(len(self.north)):
#            ax.plot(self.taps[self.north[i]].y, self.taps[self.north[i]].z, 'k+', markersize=markersize)
#            ax.axis('off')
#            border = patches.Rectangle((-self.width/2.0,0.0 ),self.width , self.height, linewidth=linewidth,edgecolor='r',facecolor='none')
#            ax.add_patch(border)
#            ax.set_title('North')
#            ax.set_ylim(0.00, self.height)            
##            ax.text(self.taps[self.north[i]].y, self.taps[self.north[i]].z, self.taps[self.north[i]].tag, fontsize=7,rotation=45)
#
#            
#        ax = fig.add_subplot(1, n_plots, 2)
#        xx, yy, zz = self.__interpolate_tap_data(ngridx, ngridy, 'West')
#        ax.contour(xx, yy, zz, 15, linewidths=0.5, colors='k')
#        ctrplt = ax.contourf(xx, yy, zz, 15, cmap="RdBu_r")
##        plt.colorbar(ctrplt, boundaries=np.linspace(0, 0.3, 15))
#        for i in range(len(self.west)):
#            ax.plot(self.taps[self.west[i]].x, self.taps[self.west[i]].z, 'k+', markersize=markersize)
#            ax.axis('off')
#            border = patches.Rectangle((-self.depth/2.0,0.0 ), self.depth , self.height, linewidth=linewidth,edgecolor='r',facecolor='none')
#            ax.add_patch(border)
#            ax.set_title('West')
#            ax.set_ylim(0.00,0.75)            
##            ax.text(self.taps[self.west[i]].x, self.taps[self.west[i]].z, self.taps[self.west[i]].tag, fontsize=7,rotation=45)
#
#            
#        ax = fig.add_subplot(1, n_plots, 3)
#        xx, yy, zz = self.__interpolate_tap_data(ngridx, ngridy, 'South')
#        ax.contour(xx, yy, zz, 15, linewidths=0.5, colors='k')
#        ctrplt = ax.contourf(xx, yy, zz, 15, cmap="RdBu_r")
##        plt.colorbar(ctrplt, boundaries=np.linspace(0, 0.3, 15))
#        for i in range(len(self.south)):
#            ax.plot(self.taps[self.south[i]].y, self.taps[self.south[i]].z, 'k+', markersize=markersize)
#            ax.axis('off')
#            border = patches.Rectangle((-self.width/2.0,0.0 ), self.width , self.height, linewidth=linewidth,edgecolor='r',facecolor='none')
#            ax.add_patch(border)
#            ax.set_title('South')            
#            ax.set_ylim(0.00,0.75)            
##            ax.text(self.taps[self.south[i]].y, self.taps[self.south[i]].z, self.taps[self.south[i]].tag, fontsize=7,rotation=45)            
#
#        ax = fig.add_subplot(1, n_plots, 4)
#        xx, yy, zz = self.__interpolate_tap_data(ngridx, ngridy, 'East')
#        ax.contour(xx, yy, zz, 15, linewidths=0.5, colors='k')
#        ctrplt = ax.contourf(xx, yy, zz, 15, cmap="RdBu_r")
##        plt.colorbar(ctrplt, boundaries=np.linspace(0, 0.3, 15))
#        for i in range(len(self.east)):
#            ax.plot(self.taps[self.east[i]].x, self.taps[self.east[i]].z, 'k+', markersize=markersize)
#            ax.axis('off')
#            border = patches.Rectangle((-self.depth/2.0,0.0 ), self.depth , self.height, linewidth=linewidth,edgecolor='r',facecolor='none')
#            ax.add_patch(border)
#            ax.set_title('East')            
#            ax.set_ylim(0.00,0.75)            
##            ax.text(self.taps[self.east[i]].x, self.taps[self.east[i]].z, self.taps[self.east[i]].tag, fontsize=7,rotation=45)
#
#
##        ax = fig.add_subplot(1, n_plots, 5)
##        for i in range(len(self.top)):
##            ax.plot(self.taps[self.top[i]].x, self.taps[self.top[i]].y, 'k+', markersize=markersize)
##            ax.axis('off')
##            border = patches.Rectangle((-self.depth/2.0, -self.width/2.0), self.depth , self.width, linewidth=linewidth, edgecolor='r', facecolor='none')
##            ax.add_patch(border)
##            ax.text(self.taps[self.top[i]].x, self.taps[self.top[i]].y, self.taps[self.top[i]].tag, fontsize=7,rotation=90)
#        
#        fig.set_size_inches(30/2.54, 75/2.54)    
#        plt.tight_layout()
#        plt.show()