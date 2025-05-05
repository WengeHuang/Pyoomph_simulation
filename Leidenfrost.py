# Import pyoomph and all used equations and utilities
from pyoomph import * # basic stuff

from pyoomph.equations.multi_component import * #multi-component flow equations
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
        
        #self.add_plot("droplet/velocity",mode="arrows",transform=[None,"mirror_x"])
        #self.add_plot("air/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("droplet/velocity",mode="arrows",transform="mirror_x")
        self.add_plot("air/velocity",mode="arrows",transform="mirror_x")
        
        self.add_plot("droplet/droplet_interface",transform=[None,"mirror_x"])
        #self.add_plot("liquid/liquid_substrate",transform=[None,"mirror_x"])
        #self.add_plot("bubble/bubble_substrate",transform=[None,"mirror_x"])

        arrkey=self.add_arrow_key("bottom right",title="water mass transfer [g/m²/s]",factor=1000)

        arrkey.ymargin=0.175
        arrkey.xmargin=0.15
        #self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform=[None,"mirror_x"])
        self.add_plot("droplet/droplet_interface/masstrans_water",arrowkey=arrkey,transform="mirror_x")




        sb=self.add_scale_bar("bottom center")
        tl=self.add_time_label("bottom center")
        tl.unit="ms"
        tl.yshift+=0.05


# Mesh class. An axisymmetric bubble with liquid around on a substrate
#class SingleBubbleAxisymmMesh(GmshTemplate):
# Create a mesh
class LeidenfrostAxisymmMesh(GmshTemplate):
    def __init__(self, droplet_height: float = 1.05*milli*meter, droplet_radius: float = 1*milli*meter, channel_height: float = 10*milli*meter, channel_width: float = 10*milli*meter):
        super().__init__()
        self.channel_height=channel_height
        self.droplet_height=droplet_height
        self.droplet_radius=droplet_radius
        self.channel_width=channel_width

    def define_geometry(self):
        self.mesh_mode="tris"
        self.default_resolution=2
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
        # attach a remesher for meshing
        self.remesher = Remesher2d(self)


# Problem class. Combining the equations, meshes, etc to a full problem
class LeidenfrostProblem(Problem):
    def __init__(self):
        super().__init__()

        # Default parameters
        self.channel_height=10*milli*meter
        self.droplet_height=1.1*milli*meter
        self.droplet_radius=1*milli*meter
        self.channel_width=10*milli*meter

        self.typical_frequency=1*kilo*hertz
        
        self.sticking_coefficient=1
        
        self.Troom=20*celsius
        self.Tdroplet0=100*celsius
        self.Tsubstrate=200*celsius
        
        self.gas=Mixture(get_pure_gas("air")+0.*get_pure_gas("water"))
        #self.gas=get_pure_gas("air")
        self.liquid=get_pure_liquid("water")
        self.substrate=get_pure_solid("borosilicate")
        
        
        self.interface=None
        
        self.g=9.81*meter/second**2
        
        self.remesh_options=RemeshingOptions(max_expansion=1.2,min_expansion=0.5,min_quality_decrease=0.5)
        self.remesh_options.active=True  

        #self.alpha = self.define_global_parameter(alpha = 0.1) 

    def get_interface_properties(self):
        if self.interface is None:
            self.interface=self.liquid | self.gas
            mass_transfer_model=HertzKnudsenSchrageMassTransferModel(self.liquid,self.gas)
            mass_transfer_model.projection_space="C2"
            mass_transfer_model.sticking_coefficient=self.sticking_coefficient
            mass_transfer_model.with_prosperetti_term=True
            self.interface.set_mass_transfer_model(mass_transfer_model)
        return self.interface


    def define_problem(self):
        # Switch to an axisymmetric coordinate system
        self.set_coordinate_system("axisymmetric")
        # Select some reasonable scales
        self.set_scaling(spatial=self.droplet_radius,temporal=1*milli*second,velocity=1*meter/second)
        # And set the remaining scales by the liquid properties
        self.liquid.set_reference_scaling_to_problem(self,temperature=20*celsius)

        # The pressure scale is given by the Laplace pressure, so we must evaluate the surface tension at some constant temperature
        self.get_interface_properties()                    
        sigma0=self.liquid.evaluate_at_condition(self.interface.surface_tension,temperature=100*celsius)
        self.set_scaling(pressure=2*sigma0/self.droplet_radius)

        # the absolute pressure is the ambient pressure plus the pressure due to the Navier-Stokes
        # absolute_pressure is used e.g. in the ideal gas law for self.gas.mass_density
        self.define_named_var(absolute_pressure=1*atm + var("pressure"))

        # Add the mesh
        self+=LeidenfrostAxisymmMesh()    
        # Simplify: Use a constant viscosity (temperature and pressure independent)
        ## if we have something temperature and pressure dependent, how should we do
        self.liquid.dynamic_viscosity=self.liquid.evaluate_at_condition("dynamic_viscosity",temperature=100*celsius)  

        # Begin to assemble the equations for each of the three domains
        droplet_eqs=MeshFileOutput() # Output. PVD/VTU files can be loaded with Paraview
        air_eqs=MeshFileOutput()
        glass_eqs=MeshFileOutput()

        # Simple moving mesh dynamics
        droplet_eqs+=PseudoElasticMesh() #LaplaceSmoothedMesh()
        air_eqs+=PseudoElasticMesh()
        glass_eqs+=PseudoElasticMesh()

         # Axis of symmetry, sets e.g. velocity_x=0
        droplet_eqs+=AxisymmetryBC()@"droplet_axis"
        air_eqs+=AxisymmetryBC()@"air_axis"
        glass_eqs+=AxisymmetryBC()@"substrate_axis"




        # Add equation to the droplet
        droplet_eqs+=CompositionFlowEquations(self.liquid,isothermal=False,initial_temperature=self.Tdroplet0,with_IC=False,compo_space="C2",gravity=self.g*vector(0,-1))       
        droplet_eqs+=InitialCondition(pressure=2*sigma0/self.droplet_radius) # Start with some reasonable pressure, since the density depends on it

        # Interface: Surface tension, mass transfer, velocity connection, Stefan flow, Marangoni flow, etc.
        interf_eqs=MultiComponentNavierStokesInterface(self.interface) 
         # Temperature must be connected (but the connection must be deactivated at the contact line, will be enforced via the substrate connection)       
        interf_eqs+=ConnectFieldsAtInterface("temperature")+DirichletBC(_lagr_conn_temperature_temperature=0)@"air_axis"
        #interf_eqs+=ConnectMeshAtInterface()
        interf_eqs+=ConnectMeshAtInterface() + DirichletBC(_lagr_conn_mesh_y=0)@"air_axis" # Mesh must stay connected        
        droplet_eqs+=interf_eqs@"droplet_interface"         
        droplet_eqs+=InitialCondition(temperature=self.Tdroplet0)
        # temperature for th droplet interface
        #air_eqs+=DirichletBC(temperature=self.Tdroplet0)@"droplet_interface"

        # Navier-Stokes for the air
        air_eqs+=CompositionFlowEquations(self.gas,isothermal=False,initial_temperature=self.Troom,compo_space="C2",gravity=self.g*vector(0,-1))       
        
        #Temperature equations for the air domian
        air_eqs+=InitialCondition(temperature=self.Troom)
        air_eqs+=DirichletBC(temperature=self.Troom)@"top"
        air_eqs+=ConnectFieldsAtInterface("temperature")@"substrate_top"
        
        air_eqs+=ConnectMeshAtInterface()@"substrate_top"
        air_eqs+=NoSlipBC()@"substrate_top"

        air_eqs+=DirichletBC(velocity_x=0)@"side_wall" 

        # Substrate just a conduction equation for the temperature. This is included in the CompositionFlowEquations as well, since we set isothermal=False there
        glass_eqs+=TemperatureConductionEquation(self.substrate,space="C2")+InitialCondition(temperature=self.Tsubstrate)      
        glass_eqs+=DirichletBC(temperature=self.Tsubstrate)@"substrate_base"
        glass_eqs+=DirichletBC(mesh_y=0)@"substrate_top"
       
        # Remesh when necessary
        droplet_eqs+=RemeshWhen(self.remesh_options)
        air_eqs+=RemeshWhen(self.remesh_options)

        
        # Add all equations to the mesh domains
        self+=droplet_eqs@"droplet"+air_eqs@"air"
        self+=glass_eqs@"substrate"


    def presolve_temperature(self):
        # Start with a reasonable temperature profile by solving it before starting the simulation
        if not self.is_initialised():
            self.initialise()
        with self.select_dofs() as dofs:
            dofs.select("droplet/temperature","air/temperature","substrate/temperature","droplet/droplet_interface/masstrans_water") #,"droplet/droplet_interface/masstrans_water"
            dofs.select("droplet/droplet_interface/_lagr_conn_temperature_temperature","air/substrate_top/_lagr_conn_temperature_temperature")
            self.solve(max_newton_iterations=20)
        self.reapply_boundary_conditions()
        self.get_mesh("droplet").ignore_initial_condition=True
        self.get_mesh("air").ignore_initial_condition=True
        self.get_mesh("substrate").ignore_initial_condition=True
