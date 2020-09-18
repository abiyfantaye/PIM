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
        if self.name == 'North':
            self.plot_name = 'Front'
        if self.name == 'West':
            self.plot_name = 'Right'
        if self.name == 'South':
            self.plot_name = 'Back'
        if self.name == 'East':
            self.plot_name = 'Left'
        if self.name == 'Top':
            self.plot_name = 'Roof'

    
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
    
    def get_extent_coordinates(self):
    
        xy = np.zeros((4, 2))        

        #Point 1
        xy[0,0] = -self.width/2.0
        xy[0,1] = 0.0
        
        #Point 2
        xy[1,0] = self.width/2.0
        xy[1,1] = 0.0
        
        #Point 3
        xy[2,0] = self.width/2.0
        xy[2,1] = self.height
        
        #Point 4
        xy[3,0] = -self.width/2.0
        xy[3,1] = self.height
        
        if self.name == 'Top':
            xy[:,1] = xy[:,1] -  self.height/2.0
                                
        return xy
    
    def get_cp(self):
        n_taps  = len(self.taps)
        cp = np.zeros((n_taps,len(self.taps[0].cp)))
        
        for i in range(n_taps):
            cp[i,:] = self.taps[i].cp
                        
        return cp    
    
    def find_nearest_tap(self, tap):
        
        dist = 1.0e20
        nearest_tap = []
        
        for i in range(len(self.taps)):
            diff = self.taps[i].coord - tap.coord
            if diff.mag() < dist and tap.name != self.taps[i].name:
                dist  = diff.mag()
                nearest_tap = self.taps[i] 
        
        return nearest_tap
    
    def create_tributary_area(self):
        
        from scipy.spatial import Voronoi, voronoi_plot_2d
        from shapely.geometry import Polygon

        x,y = self.get_coordinates()
        points = np.zeros((len(x),2))
        
        for i in range(len(x)):
            points[i,0] = x[i]
            points[i,1] = y[i]
            
        #import matplotlib.pyplot as plt
        vor = Voronoi(points)
        
        #fig = voronoi_plot_2d(vor) #plots the voronoi cells 

        regions, vertices = self.voronoi_finite_polygons_2d(vor)
        
        extent = self.get_extent_coordinates()
        
        box = Polygon(extent)

        self.trib_areas = np.zeros(len(self.taps))     

        # Creat the polygons cliping the corner cells and get trib area
        for i in range(len(regions)):
            polygon = vertices[regions[i]]
            poly = Polygon(polygon)
            poly = poly.intersection(box)  # Clipping polygon
            #polygon = [p for p in poly.exterior.coords]  #           
            #plt.fill(*zip(*polygon), alpha=0.5)
            self.trib_areas[i] = poly.area
        
        #print(self.trib_areas)
        #plt.show()
    
    def face_force(self, rho, U_ref):
        cp_data = np.loadtxt(self.Cp_file_path)

        self.cp_data = np.zeros((self.tap_count, np.shape(cp_data)[1]))
        
        for i in range(self.tap_count):
            self.cp_data[i,:] = cp_data[i,:]
            
            
    def voronoi_finite_polygons_2d(self, vor, radius=None):
        """
        Reconstruct infinite voronoi regions in a 2D diagram to finite
        regions.
        Parameters
        ----------
        vor : Voronoi
            Input diagram
        radius : float, optional
            Distance to 'points at infinity'.
        Returns
        -------
        regions : list of tuples
            Indices of vertices in each revised Voronoi regions.
        vertices : list of tuples
            Coordinates for revised Voronoi vertices. Same as coordinates
            of input vertices, with 'points at infinity' appended to the
            end.
        """
                
        if vor.points.shape[1] != 2:
            raise ValueError("Requires 2D input")
    
        new_regions = []
        new_vertices = vor.vertices.tolist()
    
        center = vor.points.mean(axis=0)
        if radius is None:
            radius = vor.points.ptp().max()*2
    
        # Construct a map containing all ridges for a given point
        all_ridges = {}
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges.setdefault(p1, []).append((p2, v1, v2))
            all_ridges.setdefault(p2, []).append((p1, v1, v2))
    
        # Reconstruct infinite regions
        for p1, region in enumerate(vor.point_region):
            vertices = vor.regions[region]
    
            if all(v >= 0 for v in vertices):
                # finite region
                new_regions.append(vertices)
                continue
    
            # reconstruct a non-finite region
            ridges = all_ridges[p1]
            new_region = [v for v in vertices if v >= 0]
    
            for p2, v1, v2 in ridges:
                if v2 < 0:
                    v1, v2 = v2, v1
                if v1 >= 0:
                    # finite ridge: already in the region
                    continue
    
                # Compute the missing endpoint of an infinite ridge
    
                t = vor.points[p2] - vor.points[p1] # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal
    
                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius
    
                new_region.append(len(new_vertices))
                new_vertices.append(far_point.tolist())
    
            # sort region counterclockwise
            vs = np.asarray([new_vertices[v] for v in new_region])
            c = vs.mean(axis=0)
            angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
            new_region = np.array(new_region)[np.argsort(angles)]
    
            # finish
            new_regions.append(new_region.tolist())
    
        return new_regions, np.asarray(new_vertices)