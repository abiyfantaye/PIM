import numpy as np
from scipy import stats
import pandas as pd
from scipy import signal
from tap import Tap
from point import Point 
from face import Face 
import CWE as cwe


   
class PIM:
    """A class representing a Pressure Integration Model(PIM)"""

    def __init__(self, 
                 cp_file_name, 
                 tap_file_name,
                 wind_direction,
                 profile_file_name,
                 case_name,
                 correct_cp,
                 building_width, 
                 building_depth,                  
                 building_height, 
                 sampling_rate,
                 test_duration,
                 start_time,
                 end_time,
                 rho,
                 z0, 
                 u_ref, 
                 z_ref, 
                 gradient_height, 
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
        self.profile_file_name = profile_file_name
        self.case_name = case_name
        self.correct_cp = correct_cp
        self.building_height = building_height
        self.building_width = building_width
        self.building_depth = building_depth
        self.sampling_rate = sampling_rate
        self.test_duration = test_duration
        self.rho = rho
        self.z0 = z0
        self.u_ref = u_ref
        self.z_ref = z_ref            
        self.gradient_height = gradient_height
        self.scale = scale
        self.broken_taps = broken_taps
        self.start_time = start_time
        self.end_time = end_time
        self.faces = []
        self.taps  = []
        self.wind_profile = []
        self.ring_height = []
        self.ring_taps = []
        self._create_taps()
        self._create_faces()        
        self._fix_broken_taps() # should be called after creating the faces

    def _read_cp_data(self):
        """
        Reads the Cp data from a given text file. The cp is read in the order 
        the taps are created. 
        """
        self.time_step = 1.0/self.sampling_rate
        
        cp_raw = np.loadtxt(self.cp_file_name)
        
        if self.correct_cp:
            cp_raw *= self._get_cp_correction_factor()

        start_index = int(np.shape(cp_raw)[1]*self.start_time/self.test_duration)
        end_index = int(np.shape(cp_raw)[1]*self.end_time/self.test_duration)
        
        cp_raw = cp_raw[:,start_index:end_index]
        self.Nt = np.shape(cp_raw)[1]
        self.time = np.arange(self.start_time, self.end_time, self.time_step)

        for i in range(self.tap_count):
            self.taps[i].cp = cp_raw[i,:]
                      
    def _read_wind_profile(self):
        """
        Reads BLWTL profiles(mean velocity and turbulence intesity) for correcting 
        the Cp from the gradient height to building height. 
        """
        self.wind_profile = np.loadtxt(self.profile_file_name)       

        
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
        
        self._read_wind_profile()
        
        self._read_cp_data()
        
        #self._remove_broken_taps() #must be called before creating the faces

    def _create_faces(self):    
        """
        Creates each faces as 'North', 'South', 'East', 'West', 'Top'.
        
        """
        #create each face
        
        self.faces.append(Face('North', self.building_width,self.building_height, Point(1.0, 0.0, 0.0)))        
        self.faces.append(Face('West', self.building_depth,self.building_height, Point(0.0, 1.0, 0.0)))
        self.faces.append(Face('South', self.building_width,self.building_height, Point(-1.0, 0.0, 0.0)))
        self.faces.append(Face('East', self.building_depth,self.building_height, Point(0.0, -1.0, 0.0)))
        self.faces.append(Face('Top', self.building_width, self.building_depth, Point(0.0, 0.0, -1.0)))
        
        self.face_count  = len(self.faces)

        #Assign each tap to corresponding faces
        for i in range(self.tap_count):
            for j in range(self.face_count):
                if self.taps[i].face == self.faces[j].name:
                    self.faces[j].taps.append(self.taps[i])


    def _get_cp_correction_factor(self):
        """
        Calculates the correction factor for the cp using the square of the 
        ratio of the mean velocity at gradieint height to the building height. 
        """
        #Correct cp with velocity ratio at roof-height and gradient height 
        u_grdient = self.get_Uav(self.gradient_height)
        u_h = self.get_Uav(self.z_ref)
                
        corr = (u_grdient/u_h)**2.0
        
        return corr
    
    def _calculate_forces(self):
        """
        Calculate the forces acting on the building in x, y and z directions. 
        """
        self.forces = np.zeros((3, self.Nt))

        for face in self.faces:
            for i in range(len(face.taps)):
                force = 0.5*self.rho*(self.u_ref**2.0)*face.taps[i].cp*face.trib_areas[i]
                self.forces[0,:] += face.normal.x*force
                self.forces[1,:] += face.normal.y*force
                self.forces[2,:] += face.normal.z*force

    
    def _calculate_moments(self):
        """
        Calculate the moments acting on the building in x, y and z directions.
        The moment is calculated by the cross product of the potision vector for 
        each tap by the force calculated for the tap.
        """
        self.moments = np.zeros((3, self.Nt))

        for face in self.faces:
            for i in range(len(face.taps)):
                force = 0.5*self.rho*(self.u_ref**2.0)*face.taps[i].cp*face.trib_areas[i]
                r_x_n  = np.cross(face.taps[i].coord.array(), face.normal.array())
                self.moments[0,:] += r_x_n[0]*force
                self.moments[1,:] += r_x_n[1]*force
                self.moments[2,:] += r_x_n[2]*force
                
    def calculate_all(self):
        """
        Calls all the functions that operate on the data.
        """
        self._create_rings()
        
        for face in self.faces:
            face.create_tributary_area()
        
        self._calculate_forces()
        self._calculate_moments()
    
    
    def get_two_third_H_cp(self):
        """
        Calculates the cp at the 2/3 height of the building for comparison 
        with litrature. 
        
        Returns
        -------
        x : distance on the surface from the left corner of the North face
        Cp: Cp values 
        
        """
        
        H23 = self.building_height*(2.0/3.0)
        
        X = [] 
        CP = []
        x_old = 0.0
            
        for i in range(self.face_count-2):
            x, cp = self.faces[i].find_horizontal_taps(H23)
            x = x 
            
            for j in range(len(x)):
                X.append(x[j] + x_old)
                CP.append(cp[j,:])
            x_old += self.faces[i].width
        
        return  np.asarray(X), np.asarray(CP)
    
    def _create_rings(self):
        """
        Finds and creates rings of pressure taps.
        """
        tap_z = np.zeros(self.tap_count)
        
        for i in range(self.tap_count):
            tap_z[i] = self.taps[i].coord.z 
        
        self.ring_height = np.unique(np.sort(tap_z))
        
        for i in range(len(self.ring_height)):
            ring_taps  = []
            for j in range(self.tap_count):
                if self.ring_height[i] == self.taps[j].coord.z:
                    ring_taps.append(self.taps[j])
            self.ring_taps.append(ring_taps)
    
    def get_Uav(self, z):
        return np.interp(z, self.wind_profile[:,0], self.wind_profile[:,1])
    
    def get_Iu(self, z):
        return np.interp(z, self.wind_profile[:,0], self.wind_profile[:,2])
    
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