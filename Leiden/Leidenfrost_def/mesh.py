# Leidenfrost droplet mesh
from pyoomph import *
from pyoomph.typings import List, Optional, Union
from pyoomph.expressions.units import *
from pyoomph.meshes.remesher import *

class LeidenfrostAxisymmMesh(GmshTemplate):
    #def __init__(self, droplet_height: float = 1.1*milli*meter, droplet_radius: float = 1*milli*meter, channel_height: float = 10*milli*meter, channel_width: float = 10*milli*meter):
    def __init__(self):
        super(LeidenfrostAxisymmMesh, self).__init__()
        # attach a remesher for meshing
        #self.remesher = Remesher2d(self)
        self.remesher=LeidenfrostRemesher(self)

    def define_geometry(self):
        pr = cast(LeidenfrostProblem,self.get_problem())
        self.mesh_mode="tris"
        self.default_resolution=2
        self.channel_height=pr.channel_height
        self.droplet_height=pr.droplet_height
        self.droplet_radius=pr.droplet_radius
        self.channel_width=pr.channel_width
        # Create points
        p00, pW0= self.point(0,0, size=0.005*self.default_resolution), self.point(self.channel_width,0)
        p0H2, pWH2=self.point(0,self.channel_height,size=0.1*self.default_resolution), self.point(self.channel_width,self.channel_height,size=2*self.default_resolution)
        p0DH, pRDHD,p0DHplus,p0DHminus=self.point(0,self.droplet_height),self.point(self.droplet_radius,self.droplet_height, size=0.025*self.default_resolution),self.point(0,self.droplet_height+self.droplet_radius, size=0.025*self.default_resolution),self.point(0,self.droplet_height-self.droplet_radius, size=0.005*self.default_resolution)
        p0S, pWS = self.point(0, -0.1*self.channel_height, size=0.1*self.default_resolution), self.point(self.channel_width,-0.1*self.channel_height)# ,size=0.025*self.default_resolution
        # Create lines
        self.create_lines(p0DHminus,"air_axis",p00,"substrate_top",pW0,"side_wall",pWH2,"top",p0H2,"air_axis",p0DHplus)
        self.create_lines(p0DHminus,"droplet_axis",p0DHplus)
        self.circle_arc(pRDHD,p0DHplus,center=p0DH,name="droplet_interface")
        self.circle_arc(pRDHD,p0DHminus,center=p0DH,name="droplet_interface")
        self.create_lines(p00,"substrate_axis",p0S,"substrate_base",pWS, "substrate_wall", pW0)
        # Create surfaces
        self.plane_surface("droplet_axis","droplet_interface",name="droplet")
        self.plane_surface("air_axis","top","side_wall","substrate_top","droplet_interface",name="air")
        self.plane_surface("substrate_axis","substrate_base","substrate_wall","substrate_top",name="substrate")

