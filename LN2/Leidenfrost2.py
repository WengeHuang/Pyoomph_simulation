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
class LeidenfrostPlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr = cast(LeidenfrostProblem, self.get_problem())
        xrange = 3*pr.droplet_radius # view range
        self.background_color = "darkgrey"
        self.set_view(-xrange, -0.2*xrange, xrange, xrange)  # -x_min, -ymin, x_max, y_max of the view window
        cb_T = self.add_colorbar("temperature [°C]", offset=-273.15, position = "bottom left")
        cb_u = self.add_colorbar("velocity [m/s]", position = "bottom right", cmap = "viridis")

        #showing mesh
        self.add_plot("droplet")
        self.add_plot("air")
        self.add_plot("substrate")


        self.add_plot("droplet/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("air/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("substrate/temperature",colorbar=cb_T,transform="mirror_x")
        
        self.add_plot("droplet/velocity",colorbar=cb_u)
        self.add_plot("air/velocity",colorbar=cb_u)
        
        self.add_plot("droplet/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("air/velocity",mode="arrows",transform=[None,"mirror_x"])
        #self.add_plot("bubble/velocity",mode="arrows",transform="mirror_x")
        #self.add_plot("liquid/velocity",mode="arrows",transform="mirror_x")
        
        self.add_plot("droplet/droplet_interface",transform=[None,"mirror_x"])
        #self.add_plot("liquid/liquid_substrate",transform=[None,"mirror_x"])
        #self.add_plot("bubble/bubble_substrate",transform=[None,"mirror_x"])

        #arrkey=self.add_arrow_key("bottom right",title="water mass transfer [kg/m²/s]")
        #arrkey.ymargin=0.175
        #arrkey.xmargin=0.15
        #self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform=[None,"mirror_x"])
        #self.add_plot("droplet/droplet_interface/masstrans_water",arrowkey=arrkey,transform="mirror_x")
        




        sb=self.add_scale_bar("bottom center")
        tl=self.add_time_label("bottom center")
        tl.unit="us"
        tl.yshift+=0.05


# Mesh class. An axisymmetric bubble with liquid around on a substrate
#class SingleBubbleAxisymmMesh(GmshTemplate):
# Create a mesh
class LeidenfrostAxisymmMesh(GmshTemplate):
    def __init__(self, droplet_height: float = 1.5*milli*meter, droplet_radius: float = 1*milli*meter, channel_height: float = 10*milli*meter, channel_width: float = 10*milli*meter):
        super().__init__()
        self.channel_height=channel_height
        self.droplet_height=droplet_height
        self.droplet_radius=droplet_radius
        self.channel_width=channel_width

    def define_geometry(self):
        self.mesh_mode="tris"
        self.default_resolution=1
        # Create points
        p00, pW0= self.point(0,0, size=0.025*self.default_resolution), self.point(self.channel_width,0)
        p0H2, pWH2=self.point(0,self.channel_height,size=0.025*self.default_resolution), self.point(self.channel_width,self.channel_height,size=2*self.default_resolution)
        p0DH, pRDHD,p0DHplus,p0DHminus=self.point(0,self.droplet_height),self.point(self.droplet_radius,self.droplet_height, size=0.025*self.default_resolution),self.point(0,self.droplet_height+self.droplet_radius, size=0.025*self.default_resolution),self.point(0,self.droplet_height-self.droplet_radius, size=0.025*self.default_resolution)
        p0S, pWS = self.point(0, -0.1*self.channel_height), self.point(self.channel_width,-0.1*self.channel_height)# ,size=0.025*self.default_resolution

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
        # attach a remesher for meshing
        self.remesher = Remesher2d(self)

# Problem class. Combining the equations, meshes, etc to a full problem
class LeidenfrostProblem(Problem):
    def __init__(self):
        super().__init__()

        # Default parameters
        self.channel_height=10*milli*meter
        self.droplet_height=1.5*milli*meter
        self.droplet_radius=1*milli*meter
        self.channel_width=10*milli*meter

        # Define physical parameters of gas, liquid and substrate
        self.gravity=9.81*meter/second**2 * vector(0,-1) # gravity direction and strength
        self.rho_liq=1000*kilogram/meter**3  # density of liquid (temperature dependency is added later)     
        self.mu_liq=1*milli*pascal*second # dynamic viscosity of liquid       
        self.Troom=20*celsius
        self.Tbubble=100*celsius
        self.Tsubstrate=200*celsius
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