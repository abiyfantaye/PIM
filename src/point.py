"""
Represents a point in three dimention, with x, y, z coordinates.

"""
class Point:
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z       
        
    def __str__(self):
        return '(' + self.x + ',' + self.y + ',' + self.z + ')' 
    
    def __add__(self, o):
        return Point(self.x + o.x, self.y + o.y, self.z + o.z)
    
    def __sub__(self, o):
        return Point(self.x - o.x, self.y - o.y, self.z - o.z)