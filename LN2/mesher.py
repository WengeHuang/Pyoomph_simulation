from pyoomph import *
from pyoomph.typings import List, Optional, Union
from pyoomph.expressions.units import *
from pyoomph.meshes.remesher import *

# Create a mesh
class Mesh(GmshTemplate):
    def __init__(self, droplet_height: float = 1.5*milli*meter, droplet_radius: float = 1*milli*meter, channel_height: float = 10*milli*meter, channel_width: float = 10*milli*meter):
        super().__init__()
        self.channel_height=channel_height
        self.droplet_height=droplet_height
        self.droplet_radius=droplet_radius
        self.channel_width=channel_width

    def define_geometry(self):
        
        # Create points
        p00, pW0= self.point(0,0, size=0.025*self.default_resolution), self.point(self.channel_width,0)
        p0H2, pWH2=self.point(0,self.channel_height,size=0.025*self.default_resolution), self.point(self.channel_width,self.channel_height,size=2*self.default_resolution)
        p0DH, pRDHD,p0DHplus,p0DHminus=self.point(0,self.droplet_height),self.point(self.droplet_radius,self.droplet_height, size=0.025*self.default_resolution),self.point(0,self.droplet_height+self.droplet_radius, size=0.025*self.default_resolution),self.point(0,self.droplet_height-self.droplet_radius, size=0.025*self.default_resolution)

        # Create lines
        self.create_lines(p0DHminus,"air_axis",p00,"substrate",pW0,"side_wall",pWH2,"top",p0H2,"air_axis",p0DHplus)
        self.circle_arc(pRDHD,p0DHplus,center=p0DH,name="droplet")
        self.circle_arc(pRDHD,p0DHminus,center=p0DH,name="droplet")

        # Create surfaces
        self.plane_surface("air_axis","top","wall","substrate","droplet",name="air")

