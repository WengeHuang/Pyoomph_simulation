from pyoomph import *
from pyoomph.equations.navier_stokes import * # Navier-Stokes for the flow
from pyoomph.equations.ALE import * # Moving mesh equations
from pyoomph.meshes.remesher import Remesher2d # Remeshing
from pyoomph.utils.num_text_out import NumericalTextOutputFile # Tool to write a file with numbers

from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools

 
# Plotter calss. Will generate figures
class HangingDropletPlotter(MatplotlibPlotter):
    def define_plot(self):

        self.background_color = "darkgrey"
        self.set_view(-2, -2, 2, 0.5)  # -x_min, -ymin, x_max, y_max of the view window
        #self.set_view(-0.5*xrange, -0.05*xrange, 0, 0.4*xrange)
        cb_p = self.add_colorbar("temperature [°C]", position = "top right")
        cb_u = self.add_colorbar("J", position = "top left", cmap = "viridis")
 
        self.add_plot("droplet/pressure",colorbar=cb_p,transform=[None,"mirror_x"])
        self.add_plot("droplet/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("droplet/velocity",colorbar=cb_u,transform="mirror_x")
 
        #self.add_plot("liquid/J",mode="streamlines",transform="mirror_x")
 
        tl=self.add_time_label("bottom center")




# Make the mesh of the hanging hemispherical droplet with radius 1
class HangingDropletMesh(GmshTemplate):
    def define_geometry(self):
        self.default_resolution=0.05
        self.mesh_mode="tris"
        self.create_lines((0,-1),"axis",(0,0),"wall",(1,0))
        self.circle_arc((0,-1),(1,0),center=(0,0),name="interface")
        self.plane_surface("axis","wall","interface",name="droplet")
        self.remesher=Remesher2d(self) # attach a remesher

class HangingDropletProblem(Problem):
    def __init__(self):
        super(HangingDropletProblem, self).__init__()
        # Initially no gravity
        self.Bo=self.define_global_parameter(Bo=0)
        # Initial volume of a hemispherical droplet with radius R=1
        self.V=self.define_global_parameter(V=2*pi/3)

    def define_problem(self):
        self.set_coordinate_system("axisymmetric")
        # add the mesh
        self+=HangingDropletMesh()
        self+=HangingDropletPlotter()

        # Bulk equation: Navier-Stokes + Moving Mesh
        eqs = MeshFileOutput()  # Paraview output
        # Mass density can be set arbitrarly with the same results. But must be >0 to get a du/dt -> mass matrix
        eqs+=NavierStokesEquations(dynamic_viscosity=1,mass_density=1,bulkforce=self.Bo*vector(0,-1))
        eqs+=LaplaceSmoothedMesh()       

        # Boundary conditions:
        eqs+=AxisymmetryBC()@"axis"
        eqs+=(DirichletBC(mesh_y=0,mesh_x=True)+NoSlipBC())@"wall" # static no-slip wall
        eqs+=NavierStokesFreeSurface(surface_tension=1)@"interface" # free surface

        # Remeshing
        eqs+=RemeshWhen(RemeshingOptions())   

        self+=eqs@"droplet"
with HangingDropletProblem() as problem:
     # Trivial, but long way: integrate in time, extract response manually from the output
     problem.run(10,outstep=1, temporal_error=1)           