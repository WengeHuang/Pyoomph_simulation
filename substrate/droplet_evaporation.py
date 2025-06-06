#  @author Christian Diddens <c.diddens@utwente.nl>
#  @author Duarte Rocha <d.rocha@utwente.nl>
#  
#  @section LICENSE
# 
#  pyoomph - a multi-physics finite element framework based on oomph-lib and GiNaC 
#  Copyright (C) 2021-2025  Christian Diddens & Duarte Rocha
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
#  The authors may be contacted at c.diddens@utwente.nl and d.rocha@utwente.nl
#
# ========================================================================


from pyoomph import *
from pyoomph.equations.ALE import * # Moving mesh
from pyoomph.equations.navier_stokes import * # Flow
from pyoomph.expressions.units import * # units
from pyoomph.utils.dropgeom import  * # utils to calculate droplet contact angle from height and radius, etc.
from pyoomph.equations.advection_diffusion import * # for the gas diffusion
from pyoomph.meshes.remesher import * # to remesh at large distortions
from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools

# Specialize the MatplotlibPlotter for the droplet problem
class DropPlotter(MatplotlibPlotter):
    # This method defines the entire plot
    def define_plot(self):    
        p=self.get_problem() # we can get access to the problem by get_problem()
        r=p.droplet_radius # and hence can access the droplet base radius
        xrange=1.5*r # in x-direction, the image will be from [-1.5*r : 1.5*r ]
        self.background_color="darkgrey" # background color

        # First step: Set the field of view (xmin,ymin,xmax,ymax)
        self.set_view(-xrange,-0.165*xrange,xrange,0.9*xrange)

        # Second step: add colorbars with different colormaps at different positions
        # with 'factor', you can control the multiplicative factor (i.e. mm/s requires a factor of 1000 to convert the m/s to mm/s)
        cb_v=self.add_colorbar("velocity [mm/s]",cmap="seismic",position="bottom right",factor=1000)
        cb_vap=self.add_colorbar("water vapor [g/m^3]",cmap="Blues",position="top right",factor=1000)

        # Now, we can add all kinds of plots
        # plot the velocity (magnitude, since it is a vector) of the droplet domain (on both sides)
        self.add_plot("droplet/velocity",colorbar=cb_v,transform=["mirror_x",None])
        # add velocity arrows on both sides
        self.add_plot("droplet/velocity", mode="arrows",linecolor="green",transform=["mirror_x",None])
        
        # Plot the vapor in the gas phase
        self.add_plot("gas/c_vap",colorbar=cb_vap,transform=["mirror_x",None])

        # at the interface lines
        self.add_plot("droplet/droplet_gas",linecolor="yellow",transform=["mirror_x",None])
        self.add_plot("droplet/droplet_substrate",transform=["mirror_x",None])
        self.add_plot("gas/gas_substrate",transform=["mirror_x",None])

        # For the evaporation rate, we require an arrow key, again with a factor 1000, since we convert from kg to g per m^2s
        ak_evap=self.add_arrow_key(position="top center",title="evap. rate water [g/(m^2s)]",factor=1000)
        # We can hide instances by setting invisible=True:
        #       ak_evap.invisible=True
        # or we can move it a bit relative to the "position" by the xmargin and ymargin
        # or with xshift and yshift. All are in graph coordinates, i.e. 1 means the width/height of the entire image
        #ak_evap.ymargin+=0.2
        #ak_evap.xmargin *= 2

        #ak_evap.rangemode="fixed"
        #ak_evap.set_range(0.1)


        ak_evap.ymargin=0.175
        ak_evap.xmargin=0.15

        # add the evaporation arrows at the interface, both sides
        arrs=self.add_plot("droplet/droplet_gas/evap_rate",arrowkey=ak_evap,transform=["mirror_x",None])

        # and a time label and a scale bar
        self.add_time_label(position="top left")
        self.add_scale_bar(position="bottom center").textsize*=0.7

# A mesh of an axisymmetric droplet surrounded by gas
class DropletWithGasMesh(GmshTemplate):
    def __init__(self,droplet_radius,droplet_height,gas_radius,cl_resolution_factor=0.1,gas_resolution_factor=50):
        super(DropletWithGasMesh, self).__init__()
        self.droplet_radius,self.droplet_height,self.gas_radius=droplet_radius,droplet_height,gas_radius
        self.cl_resolution_factor,self.gas_resolution_factor=cl_resolution_factor,gas_resolution_factor
        self.default_resolution=0.025
        #self.mesh_mode="tris"

    def define_geometry(self):
        # finer resolution at the contact line, lower at the gas far field
        res_cl=self.cl_resolution_factor*self.default_resolution
        res_gas=self.gas_resolution_factor*self.default_resolution
        # droplet
        p00=self.point(0,0) # origin
        pr0 = self.point(self.droplet_radius, 0,size=res_cl) # contact line
        p0h=self.point(0,self.droplet_height) # zenith
        self.circle_arc(pr0,p0h,through_point=self.point(-self.droplet_radius,0),name="droplet_gas") # curved interface
        self.create_lines(p0h,"droplet_axisymm",p00,"droplet_substrate",pr0)
        self.plane_surface("droplet_gas","droplet_axisymm","droplet_substrate",name="droplet") # droplet domain
        # gas dome
        pR0=self.point(self.gas_radius,0,size=res_gas)
        p0R = self.point(0,self.gas_radius,size=res_gas)
        self.circle_arc(pR0,p0R,center=p00,name="gas_infinity")
        self.line(p0h,p0R,name="gas_axisymm")
        self.line(pr0,pR0,name="gas_substrate")
        self.plane_surface("gas_substrate","gas_infinity","gas_axisymm","droplet_gas",name="gas") # gas domain


class EvaporatingDroplet(Problem):
    def __init__(self):
        super(EvaporatingDroplet, self).__init__()
        # Droplet properties
        self.droplet_radius=0.5*milli*meter # base radius
        self.droplet_height=0.2*milli*meter # apex height
        self.droplet_density=1000*kilogram/meter**3 # mass density
        self.droplet_viscosity=1*milli*pascal*second # dyn. viscosity
        self.sliplength = 1 * micro * meter # slip length

        # Gas properties
        self.gas_radius=5*milli*meter
        self.vapor_diffusivity=25.55e-6*meter**2/second
        self.c_sat=17.3*gram/meter**3 # saturated partial density of vapor
        self.c_infty=0.5*self.c_sat # partial vapor density far away

        # Interface and contact line properties
        self.surface_tension=72*milli*newton/meter # surface tension
        self.pinned_contact_line=True
        self.contact_angle=None # will be calculated from the height and radius if not set

        # Bind the evaporation rate
        c_vap = var("c_vap", domain="gas")
        n = var("normal")
        self.evap_rate = -self.vapor_diffusivity * dot(grad(c_vap), n)


    def define_problem(self):
        # Settings: Axisymmetric and typical scales
        self.set_coordinate_system("axisymmetric")
        self.set_scaling(temporal=1*second,spatial=self.droplet_radius)
        self.set_scaling(velocity=scale_factor("spatial")/scale_factor("temporal"))
        self.set_scaling(pressure=self.surface_tension/scale_factor("spatial"))
        self.set_scaling(c_vap=self.c_sat)

        # Add the mesh
        mesh=DropletWithGasMesh(self.droplet_radius,self.droplet_height,self.gas_radius)
        mesh.remesher=Remesher2d(mesh) # add remeshing possibility
        self.add_mesh(mesh)

        # Calculate the contact angle if not set
        if self.contact_angle is None:
            self.contact_angle= DropletGeometry(base_radius=self.droplet_radius, apex_height=self.droplet_height).contact_angle

        # Droplet equations
        d_eqs=MeshFileOutput() # Output
        d_eqs+=PseudoElasticMesh() # Mesh motion
        d_eqs+=NavierStokesEquations(mass_density=self.droplet_density,dynamic_viscosity=self.droplet_viscosity) # flow
        d_eqs+=DirichletBC(mesh_x=0,velocity_x=0)@"droplet_axisymm" # symmetry axis
        d_eqs += DirichletBC(mesh_y=0, velocity_y=0) @ "droplet_substrate" # allow slip, but fix mesh
        d_eqs += NavierStokesSlipLength(self.sliplength)@ "droplet_substrate" # limit slip by slip length
        #d_eqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension,mass_transfer_rate=self.evap_rate)@"droplet_gas" # Free surface equation
        d_eqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension)@"droplet_gas"
        d_eqs += ConnectMeshAtInterface() @ "droplet_gas" # connect the gas mesh to co-move

        # Different contact line dynamics
        if self.pinned_contact_line: # if pinned
            # Pinned contact line means mesh_x is fixed.
            # We enforce partial_t(mesh_x,ALE=False)=0 by adjusting the radial velocity at the contact line
            cl_constraint=mesh_velocity()[0]-0
            d_eqs+=EnforcedBC(velocity_x=cl_constraint)@"droplet_gas/droplet_substrate"
        else:
            d_eqs += NavierStokesContactAngle(contact_angle=self.contact_angle) @ "droplet_gas/droplet_substrate"  # and constant contact angle

        # Gas equations
        g_eqs=MeshFileOutput() # output
        g_eqs+=PseudoElasticMesh() # mesh motion
        g_eqs+=AdvectionDiffusionEquations(fieldnames="c_vap",diffusivity=self.vapor_diffusivity,wind=0) # diffusion equation
        g_eqs += InitialCondition(c_vap=self.c_infty)
        g_eqs+=DirichletBC(mesh_x=0)@"gas_axisymm" # fixed mesh coordinates at the boundaries
        g_eqs+=DirichletBC(mesh_y=0)@"gas_substrate"
        g_eqs+=DirichletBC(mesh_x=True,mesh_y=True)@"gas_infinity"
        g_eqs+=DirichletBC(c_vap=self.c_sat)@"droplet_gas"
        g_eqs+=AdvectionDiffusionInfinity(c_vap=self.c_infty)@"gas_infinity"

        # Control remeshing
        d_eqs += RemeshWhen(RemeshingOptions(max_expansion=1.5, min_expansion=0.7))
        g_eqs+=RemeshWhen(RemeshingOptions(max_expansion=1.5,min_expansion=0.7))

        # Output of the volume evolution
        d_eqs+=IntegralObservables(volume=1)
        d_eqs+=IntegralObservableOutput(filename="EVO_droplet")

        # Also output the interface data, along with the evaporation rate
        d_eqs+=(LocalExpressions(evap_rate=self.evap_rate)+MeshFileOutput())@"droplet_gas"

        self.add_equations(d_eqs@"droplet"+g_eqs@"gas")


if __name__=="__main__":
    with EvaporatingDroplet() as problem:
        problem.plotter=DropPlotter(problem)
        problem.run(100*second,outstep=10*second,temporal_error=1)