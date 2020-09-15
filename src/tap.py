"""
This class represents a tap in a pressure integration model 
"""

class Tap:
    
    def __init__(self, index, name, coord, face):
        self.index = index      #index used to refer it in global tap array
        self.name = name        #label used to name the tap in the presure measurment file
        self.coord = coord      #coordinate of the tap in global coordinate system.  
        self.face  = face       #represents the side on with the probe is located
        self.trib_area = 0.0    #tributary area asociated with tap.

    def print_tap_info(self):
        print('Name = ' + self.name)
        print('Index = %d' % self.index)
        print('x = %f' % self.coord.x)
        print('y = %f' % self.coord.y)
        print('z = %f' % self.coord.z)        
        print('Face = ' +  self.face)
