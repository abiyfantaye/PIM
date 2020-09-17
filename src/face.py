import numpy as np

"""
Represents face of the building model. 

"""
  
class Face:
    
    def __init__(self, name, width, height, normal):
        self.name = name        #name of the face
        self.width  = width     #width of the face 
        self.height = height    #height of the face
        self.normal = normal    #unit vector in wall normal direction
        self.trib_areas = []    #tributary area for each tap         
        self.cp = []            #
        self.taps = []
    
    def get_coordinates(self):
        n_taps  = len(self.taps)
        x = np.zeros(n_taps)
        y = np.zeros(n_taps)
        
        for i in range(n_taps):
            if self.name == 'North':
                x[i] = -self.taps[i].coord.y 
                y[i] = self.taps[i].coord.z 
            if self.name == 'West':
                x[i] = self.taps[i].coord.x  
                y[i] = self.taps[i].coord.z 
            if self.name == 'South':
                x[i] = self.taps[i].coord.y 
                y[i] = self.taps[i].coord.z 
            if self.name == 'East':
                x[i] = -self.taps[i].coord.x 
                y[i] = self.taps[i].coord.z 
            if self.name == 'Top':
                x[i] = self.taps[i].coord.y 
                y[i] = -self.taps[i].coord.x
                        
        return x, y
    
    def get_cp(self):
        n_taps  = len(self.taps)
        cp = np.zeros((n_taps,len(self.taps[0].cp)))
        
        for i in range(n_taps):
            cp[i,:] = self.taps[i].cp
                        
        return cp        
    
    def __calculate_tributary_area(self):
        cp_data = np.loadtxt(self.Cp_file_path)

        self.cp_data = np.zeros((self.tap_count, np.shape(cp_data)[1]))
        
        for i in range(self.tap_count):
            self.cp_data[i,:] = cp_data[i,:]
    
    def face_force(self, rho, U_ref):
        cp_data = np.loadtxt(self.Cp_file_path)

        self.cp_data = np.zeros((self.tap_count, np.shape(cp_data)[1]))
        
        for i in range(self.tap_count):
            self.cp_data[i,:] = cp_data[i,:]