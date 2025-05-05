from pyoomph import *
from pyoomph.typings import List, Optional, Union
from pyoomph.expressions.units import *
from pyoomph.meshes.remesher import *

# Create a mesh
class Mesh(GmshTemplate):
    def __init__(self, bubble_height: float = 125*micro*meter, bubble_radius: float = 10*micro*meter, channel_height: float = 300*micro*meter, channel_width: float = 250*micro*meter):
        super().__init__()
        self.channel_height=channel_height
        self.bubble_height=bubble_height
        self.bubble_radius=bubble_radius
        self.channel_width=channel_width

    def define_geometry(self):
        
        # Create points
        p00, pW0= self.point(0,0, size=0.025*self.default_resolution), self.point(self.channel_width,0)
        p0H2, pWH2=self.point(0,self.channel_height,size=0.025*self.default_resolution), self.point(self.channel_width,self.channel_height,size=2*self.default_resolution)
        p0BH, pRBHB,p0BHplus,p0BHminus=self.point(0,self.bubble_height),self.point(self.bubble_radius,self.bubble_height, size=0.025*self.default_resolution),self.point(0,self.bubble_height+self.bubble_radius, size=0.025*self.default_resolution),self.point(0,self.bubble_height-self.bubble_radius, size=0.025*self.default_resolution)

        # Create lines
        self.create_lines(p0BHminus,"liquid_axis",p00,"ice_front",pW0,"liquid_wall",pWH2,"liquid_top",p0H2,"liquid_axis",p0BHplus)
        self.circle_arc(pRBHB,p0BHplus,center=p0BH,name="bubble")
        self.circle_arc(pRBHB,p0BHminus,center=p0BH,name="bubble")

        # Create surfaces
        self.plane_surface("liquid_axis","liquid_top","liquid_wall","ice_front","bubble",name="liquid")