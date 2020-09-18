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
                 sampling_rate,
                 test_duration,
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
        self.sampling_rate = sampling_rate
        self.test_duration = test_duration
        self.z0 = z0
        self.u_ref = u_ref
        self.z_ref = z_ref            
        self.gradient_height = gradient_height
        self.gradient_wind_speed = gradient_wind_speed
        self.scale = scale
        self.broken_taps = broken_taps
        self.faces = []
        self.taps  = []
        self.cp_raw = []
        self._create_taps()
        self._create_faces()
        self._fix_broken_taps() # should be called after creating the faces
#        self.__correct_cp_to_building_height()

    def _read_cp_data(self):
        cp_raw = np.loadtxt(self.cp_file_name)
        
        self.Nt = np.shape(cp_raw)[1]
        self.time_step = 1.0/self.sampling_rate
        self.time = np.arange(0.0, self.test_duration, self.time_step)

        for i in range(self.tap_count):
            self.taps[i].cp = cp_raw[i,:]
                    
            

    def _create_taps(self):
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
        
        self._read_cp_data()
        
        #self._remove_broken_taps() #must be called before creating the faces

    def _create_faces(self):    
        """
        Creates each faces as 'North', 'South', 'East', 'West', 'Top'.
        
        """
        #create each face
        
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


#    def _correct_cp_to_building_height(self):
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
         

    def _find_tap(self, name):
                
        for tap in self.taps:
            if tap.name == name:
                return tap
            
    def _remove_broken_taps(self):
        """
        If any tap is broken, this fucntion removes it.
        """
        if self.broken_taps:
            for i in range(len(self.broken_taps)):
                tap = self._find_tap(self.broken_taps[i])
                self.taps.remove(tap)
                self.tap_count -= 1

    def _fix_broken_taps(self):
        """
        If any tap is broken, assignes it's value from a nearest tap.
        """
        if self.broken_taps:
            for i in range(len(self.broken_taps)):                    
                tap = self._find_tap(self.broken_taps[i])
                for face in self.faces:
                    if tap.face == face.name:
                        near_tap  = face.find_nearest_tap(tap)
                        tap.cp = near_tap.cp