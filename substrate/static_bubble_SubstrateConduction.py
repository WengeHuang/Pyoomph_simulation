# Import pyoomph and all used equations and utilities
from pyoomph import * # basic stuff

from pyoomph.equations.advection_diffusion import *
#from pyoomph.equations.multi_component import * #multi-component flow equations
from pyoomph.equations.ALE import * # moving mesh equations
from pyoomph.equations.navier_stokes import * #Navier-Stokes stuff, e.g. slip length 
from pyoomph.equations.contact_angle import * #Contact angle equations

from pyoomph.utils.dropgeom import DropletGeometry # droplet geometry calculations

from pyoomph.materials import * # Material API
from pyoomph.materials.mass_transfer import * # Mass transfer models
import pyoomph.materials.default_materials # Default materials like water, water vapor, glass

from pyoomph.meshes.remesher import Remesher2d # Remeshing on larger deformation
from pyoomph.meshes.zeta import * # Zeta coordinate interpolation upon remeshing

from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools
from pyoomph.equations.poisson import *

# Plotter calss. Will generate figures
class SingleBubblePlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr = cast(SingleBubbleProblem, self.get_problem())
        xrange = 3*pr.R0 # view range
        self.background_color = "darkgrey"
        self.set_view(-xrange, -0.3*xrange, xrange, xrange)  # -x_min, -ymin, x_max, y_max of the view window
        #self.set_view(-0.5*xrange, -0.05*xrange, 0, 0.4*xrange)
        cb_T = self.add_colorbar("temperature [°C]", offset=-273.15, position = "top left")
        cb_u = self.add_colorbar("velocity [m/s]", position = "bottom right", cmap = "viridis")

        #showing mesh
        self.add_plot("liquid")
        self.add_plot("bubble")
        self.add_plot("substrate")


        self.add_plot("liquid/temperature", colorbar=cb_T, transform = "mirror_x")
        #self.add_plot("bubble/c", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("substrate/temperature",colorbar=cb_T,transform="mirror_x")
        
        self.add_plot("liquid/velocity",colorbar=cb_u)
        #self.add_plot("bubble/velocity",colorbar=cb_u)
        #self.add_plot("bubble/c",colorbar=cb_u)
        
        #self.add_plot("bubble/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("liquid/velocity",mode="arrows",transform="mirror_x")
        #self.add_plot("bubble/velocity",mode="arrows",transform="mirror_x")
        #self.add_plot("liquid/velocity",mode="arrows",transform="mirror_x")
        
        #self.add_plot("liquid/interface",transform=[None,"mirror_x"])
        #self.add_plot("liquid/liquid_substrate",transform=[None,"mirror_x"])
        #self.add_plot("bubble/bubble_substrate",transform=[None,"mirror_x"])

        #arrkey=self.add_arrow_key("bottom right",title="water mass transfer [kg/m²/s]")
        #arrkey.ymargin=0.175
        #arrkey.xmargin=0.15
        #self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform=[None,"mirror_x"])
        #self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform="mirror_x")

# Mesh class. An axisymmetric bubble with liquid around on a substrate
class SingleBubbleAxisymmMesh(GmshTemplate):
    def define_geometry(self):
        pr=cast(SingleBubbleProblem,self.get_problem())
        geom=DropletGeometry(curv_radius=pr.R0,contact_angle=pr.contact_angle)
        self.mesh_mode="tris"
        self.default_resolution=0.1
        cl_factor=0.1
        p00=self.point(0,0)
        pr0=self.point(geom.base_radius,0,size=self.default_resolution*cl_factor)
        p0h=self.point(0,geom.apex_height)
        
        self.create_lines(p0h,"bubble_axis",p00,"bubble_substrate",pr0)
        self.circle_arc(pr0,p0h,through_point=(-geom.base_radius,0),name="interface")
        self.plane_surface("interface","bubble_substrate","bubble_axis",name="bubble")
        
        far_resolution=self.default_resolution*(pr.Lhost/pr.R0)
        pL0=self.point(pr.Lhost,0,size=far_resolution)
        p0L=self.point(0,pr.Lhost,size=far_resolution)
        pLL=self.point(pr.Lhost,pr.Lhost,size=far_resolution)
        self.create_lines(pr0,"liquid_substrate",pL0,"liquid_side",pLL,"liquid_top",p0L,"liquid_axis",p0h)
        self.plane_surface("liquid_axis","liquid_substrate","liquid_side","interface","liquid_top",name="liquid")

        # add a substrate with thickness
        p0S=self.point(0,-pr.thickness, size=far_resolution)
        pLS=self.point(pr.Lhost,-pr.thickness,size=far_resolution)
        self.create_lines(p00,"substrate_axis",p0S,"substrate_base",pLS,"substrate_side",pL0)
        self.plane_surface("substrate_axis","substrate_base","substrate_side","liquid_substrate","bubble_substrate",name="substrate")

        #p0S=self.point(0,-(pr.Lhost+pr.thichness),size=far_resolution)
        #self.circle_arc(p0S,pL0,center=p00,name="substrate_far")
        #self.line(p00,p0S,name="substrate_axis")
        #self.plane_surface("substrate_axis","substrate_far","liquid_substrate","bubble_substrate",name="substrate")


        # attach a remesher for meshing
        self.remesher = Remesher2d(self)

# Problem class. Combining the equations, meshes, etc to a full problem
class SingleBubbleProblem(Problem):
    def __init__(self):
        super().__init__()

        # Define dimensions
        self.R0=4*micro*meter
        self.contact_angle=160*degree      
        self.Lhost=100*micro*meter
        self.thickness=100*micro*meter

        # Define physical parameters of gas, liquid and substrate
        self.gravity=9.81*meter/second**2 * vector(0,-1) # gravity direction and strength
        self.rho_liq=1000*kilogram/meter**3  # density of liquid (temperature dependency is added later)     
        self.mu_liq=1*milli*pascal*second # dynamic viscosity of liquid       
        self.Troom=20*celsius
        self.Tsat=100*celsius
        self.Tsubstrate=120*celsius
        self.c0=0.01*kilogram/meter**3 # top boundary concentration
        self.csat=0.0286*kilogram/meter**3 # saturation concentration
        self.thermal_diffusivity = 1.4558e-7*meter**2/second

        # Bubble properties
        self.surface_tension=0.0721*newton/meter # surface tension (temperature dependency is added later)
        self.sigma0=0.0721*newton/meter # surface tension at c=0

        # Gas in bubble properties
        self.R_air=287.058*joule/(kilogram*kelvin) # gas constant of air
        self.rho_air_init=1.2754*kilogram/meter**3 # density of air
        self.mu_air=1.729e-5*pascal*second # dynamic viscosity of air
        self.diffusivity=1e-9*meter**2/second # diffusivity of air in water
        self.pressure_0=self.rho_air_init*self.R_air*self.Troom # pressure initially
        # Substrate properties
        self.molar_mass = 122.0632 * gram / mol  ##TODO: correct value
        # https://en.wikipedia.org/wiki/Borosilicate_glass
        self.mass_density = 2230.0 * kilogram / (meter**3)        
        self.specific_heat_capacity = 830 * joule / (kilogram * kelvin)
        # https://www.design1st.com/Design-Resource-Library/engineering_data/MaterialPropertiesTables.pdf
        self.thermal_conductivity = 1.13 * watt / (meter * kelvin)    
        self.thermal_diffusivity_s = self.thermal_conductivity/ (self.mass_density*self.specific_heat_capacity)

        # Contact line define
        #self.pinned_contact_line=True
        self.pinned_contact_line=False
        
        self.sliplength=0.01*micro*meter
       
        self.remesh_options=RemeshingOptions()
        self.remesh_options.active=True 

        # Values at current time step
        self.rho_air=var("gas_mass", domain="bubble")/var("gas_volume", domain="bubble") # density of air
        self.mass_transfer_rate = self.diffusivity*dot(var("normal"),grad(var("c", domain="liquid"))) # mass transfer rate


    
    def define_problem(self):
        
        # Redefine physical parameters according to the model
        TKelvin=var("temperature")/scale_factor("temperature")
        self.surface_tension+=0.07275*0.002*(TKelvin-291.0)*newton/meter # surface tension
        self.sigma0+=0.07275*0.002*(self.Tsubstrate/scale_factor("temperature")-291.0)*newton/meter # surface tension
        self.rho_liq+=(1000 * ((TKelvin + 15.7914) / (508929 * (TKelvin - 205.02037)) * (TKelvin - 277.1363) ** 2))*kilogram/meter**3 # density of liquid
    

        # Switch to an axisymmetric coordinate system
        self.set_coordinate_system("axisymmetric")

        # Scales for nondimensionalization
        self.set_scaling(spatial=self.R0)
        self.set_scaling(temporal=1*micro*second)
        self.set_scaling(velocity=1*meter/second)
        self.set_scaling(pressure=self.sigma0/scale_factor("spatial"))
        self.set_scaling(c=self.c0)
        self.set_scaling(temperature=self.Troom)

        self.define_named_var(absolute_pressure=1*atm + var("pressure"))

        # add the mesh
        self+=SingleBubbleAxisymmMesh()
        self+=SingleBubblePlotter()

        geom=DropletGeometry(curv_radius=self.R0,contact_angle=self.contact_angle)


        #Mesh and output
        # Begin to assemble the equations for each of the three domains
        liquid_eqs=MeshFileOutput()
        bubble_eqs=MeshFileOutput()
        glass_eqs=MeshFileOutput()


        # Describe velocity equations for liquid
        liquid_eqs+=NavierStokesEquations(mass_density=self.rho_liq,dynamic_viscosity=self.mu_liq,gravity=self.gravity)
        liquid_eqs+=AxisymmetryBC()@"liquid_axis"
        liquid_eqs+=LaplaceSmoothedMesh()
        liquid_eqs+=DirichletBC(velocity_x=0)@"liquid_side"
        liquid_eqs+=DirichletBC(mesh_y=0,velocity_y=0)@"liquid_substrate"
        liquid_eqs+=NoSlipBC()@"liquid_substrate"    
        liquid_eqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension, static_interface=True)@"interface"  #, static_interface=True
        #liquid_eqs+=ConnectMeshAtInterface() @ "interface" # connect the gas mesh to co-move
        #liquid_eqs+=ConnectMeshAtInterface() @ "liquid_substrate" # connect the gas mesh to co-move

        # Add temperature effects
        liquid_eqs+=AdvectionDiffusionEquations("temperature", diffusivity=self.thermal_diffusivity, space="C1", wind=var("velocity"))
        liquid_eqs+=AdvectionDiffusionInfinity(temperature=self.Troom)@"liquid_side"
        liquid_eqs+=AdvectionDiffusionInfinity(temperature=self.Troom)@"liquid_top"
        liquid_eqs+=InitialCondition(temperature=self.Troom)
        #liquid_eqs+=DirichletBC(temperature=self.Tsubstrate)@"liquid_substrate"
        liquid_eqs+=ConnectFieldsAtInterface("temperature")@"liquid_substrate"


        # leqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C1", wind=0)
        #liquid_eqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C1", wind=var("velocity"))
        #liquid_eqs+=AdvectionDiffusionInfinity(c=self.c0)@"liquid_top"
        #liquid_eqs+=AdvectionDiffusionInfinity(c=self.c0)@"liquid_side"
        #liquid_eqs+=DirichletBC(c=self.csat)@"interface"
        #liquid_eqs+=InitialCondition(c=self.c0)


        # Equations for bubble
        #bubble_eqs+=LaplaceSmoothedMesh()
        #bubble_eqs+=AxisymmetryBC()@"bubble_axis"
        #bubble_eqs+=DirichletBC(mesh_y=0)@"bubble_substrate"
        #bubble_eqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C1", wind=0)
        #bubble_eqs+=InitialCondition(c=self.c0)

        # Remesh when necessary
        #liquid_eqs+=RemeshWhen(self.remesh_options)
        #bubble_eqs+=RemeshWhen(self.remesh_options)
        #glass_eqs+=RemeshWhen(self.remesh_options)

        # Equation 
        #glass_eqs+=LaplaceSmoothedMesh()
        glass_eqs+=AdvectionDiffusionEquations("temperature",diffusivity=self.thermal_diffusivity_s, space="C1", wind=0)
        glass_eqs+=DirichletBC(temperature=self.Tsubstrate)@"substrate_base"
        glass_eqs+=AxisymmetryBC()@"substrate_axis"


        # Add equations
        # Add all equations to the mesh domains
        self+=bubble_eqs@"bubble"+liquid_eqs@"liquid"
        self+=glass_eqs@"substrate"

    def presolve_temperature(self):
        # Start with a reasonable temperature profile by solving it before starting the simulation
        if not self.is_initialised():
            self.initialise()
        with self.select_dofs() as dofs:
            dofs.select("liquid/temperature","substrate/temperature", "liquid/liquid_substrate/_lagr_conn_temperature_temperature") #,"substrate/temperature"
            #dofs.select("liquid/c")
            self.solve()
        self.reapply_boundary_conditions()
        self.get_mesh("liquid").ignore_initial_condition=True


with SingleBubbleProblem() as p:
    p.presolve_temperature()
    p.run(10*micro*second, outstep=0.1*micro*second, temporal_error=1)
