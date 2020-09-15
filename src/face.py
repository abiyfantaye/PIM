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