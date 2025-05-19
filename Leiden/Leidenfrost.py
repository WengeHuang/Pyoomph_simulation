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


from pyoomph.materials.mass_transfer import DifferenceDrivenMassTransferModel # Mass transfer models

# Plotter calss. Will generate figures
class LeidenfrostPlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr = cast(LeidenfrostProblem, self.get_problem())
        xrange = 3*pr.droplet_radius # view range
        self.background_color = "darkgrey"
        #self.set_view(-0.01*xrange, -0.02*xrange, 0.2*xrange, 0.09*xrange)  # -x_min, -ymin, x_max, y_max of the view window
        self.set_view(-xrange, -0.2*xrange, xrange, xrange)
        cb_T = self.add_colorbar("temperature [°C]", offset=-273.15, position = "bottom left",cmap="plasma")
        cb_u = self.add_colorbar("velocity [m/s]", position = "bottom right", cmap = "coolwarm") #cmap = "coolwarm" cmap = "viridis"

        #showing mesh
        #self.add_plot("droplet")
        #self.add_plot("air")
        #self.add_plot("substrate")

        self.add_plot("droplet/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("air/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("substrate/temperature",colorbar=cb_T,transform="mirror_x")
        
        self.add_plot("droplet/velocity",colorbar=cb_u)
        self.add_plot("air/velocity",colorbar=cb_u)
        
        self.add_plot("droplet/velocity",mode="arrows",transform=[None,"mirror_x"], linecolor="white")
        #self.add_plot("air/velocity",mode="arrows",transform="mirror_x")
        self.add_plot("air/velocity",mode="arrows",transform=[None,"mirror_x"])

        self.add_plot("droplet/droplet_interface",transform=[None,"mirror_x"])
        self.add_plot("substrate/substrate_top",transform=[None,"mirror_x"])
        arrkey=self.add_arrow_key("top right",title="water mass transfer [g/m²/s]")
        arrkey.ymargin=0.175
        arrkey.xmargin=0.15
        arrkey.rangemode="fixed"
        arrkey.set_range(0.1)
        self.add_plot("droplet/droplet_interface/masstrans_water",arrowkey=arrkey,transform="mirror_x")




        sb=self.add_scale_bar("bottom center")
        sb.yshift+=0.05
        tl=self.add_time_label("bottom center")
        tl.unit="ms"
        tl.yshift-=0.04


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
        #geometry inf
        self.channel_height=pr.channel_height
        self.droplet_height=pr.droplet_height
        self.droplet_radius=pr.droplet_radius
        self.channel_width=pr.channel_width
        #mesh size inf
        self.mesh_size_bottom = 0.01
        self.mesh_size_center = 0.05
        self.mesh_size_top = 0.05
        self.mesh_size_far = 0.5

        # Create points
        p00 = self.point(0,0, size=self.mesh_size_bottom*self.default_resolution)
        # Droplet center
        p0DH = self.point(0,self.droplet_height, size=self.mesh_size_center*self.default_resolution)  
        # Droplet equator    
        pRDHD = self.point(self.droplet_radius,self.droplet_height, size=self.mesh_size_top*self.default_resolution)
        # Droplet top
        p0DHplus = self.point(0,self.droplet_height+self.droplet_radius, size=self.mesh_size_top*self.default_resolution)
        # Droplet bottom
        p0DHminus = self.point(0,self.droplet_height-self.droplet_radius, size=self.mesh_size_bottom*self.default_resolution)

        # Gas side
        pW0 = self.point(self.channel_width,0, size=self.mesh_size_far*self.default_resolution)
        # Gas top
        p0H2 = self.point(0,self.channel_height,size=self.mesh_size_far*self.default_resolution)
        # Gas corner
        pWH2 = self.point(self.channel_width,self.channel_height,size=self.mesh_size_far*self.default_resolution)
        # Substrate bottom
        p0S = self.point(0, -0.1*self.channel_height, size=0.1*self.default_resolution)
        # Substrate corner
        pWS = self.point(self.channel_width,-0.1*self.channel_height, size=self.mesh_size_far*self.default_resolution)

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

class LeidenfrostRemesher(Remesher2d):
    def _define_geometry(self):
        super()._define_geometry()
        # geometry of the box size
        xmin, ymin, xmax, ymax = 0, -0.05, 0.8, 0.1     
        # resolution of the remesh
        mesh_in = 0.005
        mesh_out = 2                   
        import gmsh
        box = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(box, "VIn", mesh_in)   
        gmsh.model.mesh.field.setNumber(box, "VOut", mesh_out)   

        gmsh.model.mesh.field.setNumber(box, "XMin", xmin)
        gmsh.model.mesh.field.setNumber(box, "XMax", xmax)
        gmsh.model.mesh.field.setNumber(box, "YMin", ymin)
        gmsh.model.mesh.field.setNumber(box, "YMax", ymax)

        gmsh.model.mesh.field.setAsBackgroundMesh(box)

        # Add a distance field measusing the distance to the interface
        #distance=gmsh.model.mesh.field.add("Distance")
        #interface=[l.gmsh_line._id for l in self.get_line_entries_by_phys_name("substrate_top")]
        #gmsh.model.mesh.field.setNumbers(distance, "CurvesList", interface)
        #interface1 = [l.gmsh_line._id for l in self.get_line_entries_by_phys_name("droplet_interface")]
        #interface2 = [l.gmsh_line._id for l in self.get_line_entries_by_phys_name("substrate_top")]
        #interfaces = interface1 + interface2  
        #gmsh.model.mesh.field.setNumbers(distance, "CurvesList", interfaces)

        #rdist=0.5 # Blending distance
        #threshold = gmsh.model.mesh.field.add("Threshold")
        #gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
        #gmsh.model.mesh.field.setNumber(threshold, "LcMin", 0.005) # fine resolution closer to the interface
        #gmsh.model.mesh.field.setNumber(threshold, "LcMax", 2) # far away from the interface
        #gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.1*rdist) # Blending distances
        #gmsh.model.mesh.field.setNumber(threshold, "DistMax", 5*rdist)
        #gmsh.model.mesh.field.setAsBackgroundMesh(threshold)

# Problem class. Combining the equations, meshes, etc to a full problem
class LeidenfrostProblem(Problem):
    def __init__(self):
        super().__init__()

        # Default parameters
        self.channel_height=10*milli*meter
        self.droplet_height=1.5*milli*meter
        self.droplet_radius=1*milli*meter
        self.channel_width=10*milli*meter

        self.typical_frequency=1*kilo*hertz #1*kilo*hertz
        
        self.sticking_coefficient=1
        
        self.Troom=20*celsius
        self.Tdroplet0=100*celsius
        self.Tsubstrate=300*celsius
        
        self.gas=Mixture(get_pure_gas("air")+20*percent *get_pure_gas("water"))
        #self.gas=get_pure_gas("air")
        self.liquid=get_pure_liquid("water")
        self.substrate=get_pure_solid("borosilicate")
        
        
        self.interface=None
        
        self.g=9.81*meter/second**2
        
        self.remesh_options=RemeshingOptions()#max_expansion=1.2,min_expansion=0.5,min_quality_decrease=0.5
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
        droplet_eqs+=LaplaceSmoothedMesh() #LaplaceSmoothedMesh()PseudoElasticMesh()
        air_eqs+=LaplaceSmoothedMesh()
        glass_eqs+=LaplaceSmoothedMesh()

         # Axis of symmetry, sets e.g. velocity_x=0
        droplet_eqs+=AxisymmetryBC()@"droplet_axis"
        air_eqs+=AxisymmetryBC()@"air_axis"
        glass_eqs+=AxisymmetryBC()@"substrate_axis"




        # Add equation to the droplet
        droplet_eqs+=CompositionFlowEquations(self.liquid,isothermal=False,initial_temperature=self.Tdroplet0,with_IC=False,compo_space="C2",gravity=self.g*vector(0,-1))       
        droplet_eqs+=InitialCondition(pressure=2*sigma0/self.droplet_radius) # Start with some reasonable pressure, since the density depends on it

        # Interface: Surface tension, mass transfer, velocity connection, Stefan flow, Marangoni flow, etc.
        interf_eqs=MultiComponentNavierStokesInterface(self.interface, total_mass_loss_factor_inside=0) 
         # Temperature must be connected (but the connection must be deactivated at the contact line, will be enforced via the substrate connection)       
        interf_eqs+=ConnectFieldsAtInterface("temperature")+DirichletBC(_lagr_conn_temperature_temperature=0)@"air_axis"
        #interf_eqs+=ConnectMeshAtInterface()
        interf_eqs+=ConnectMeshAtInterface() + DirichletBC(_lagr_conn_mesh_y=0)@"air_axis" # Mesh must stay connected        
        droplet_eqs+=interf_eqs@"droplet_interface"         
        droplet_eqs+=InitialCondition(temperature=self.Tdroplet0)
        # constant volume constrain
        V0 = 4/3*pi*self.droplet_radius**3
        droplet_eqs+=EnforceVolumeByPressure(V0)
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

        # Measure the volume (integral over 1 ) 
        droplet_eqs+=IntegralObservables(volume=1)
        # Droplet center of mass height and velocity
        droplet_eqs+=IntegralObservables(_volume=1, _height = var("coordinate_y"), height = lambda _volume,_height : _height/_volume)
        droplet_eqs+=IntegralObservables(_velocity=var("velocity_y"), velocity = lambda _volume,_velocity : _velocity/_volume)
        ##droplet_eqs+=TextFileOutput()@"droplet_interface/droplet_axis" # all inf on the interface
        droplet_eqs+=IntegralObservableOutput("droplet_volume_height_velocity") # write it to bubble_evolution.txt
        
        # Measure the surface area (integral over 1 ) 
    #    droplet_eqs+=IntegralObservables(area=1)@"droplet_interface"
    #    droplet_eqs+=IntegralObservableOutput("surface_area")@"droplet_interface"
        
        # Measure the some data at points 
    #    only_lower=heaviside(-var("normal_y",domain=".."))  # Only take the bottom point, not the upper a step function
    #    droplet_eqs+=IntegralObservables(cartesian,h=only_lower*var("coordinate_y"),pressure=only_lower*var("pressure"),T=only_lower*var("temperature"))@"droplet_interface/droplet_axis"  # And potentially more output   
        # the output pressure var("pressure") is just the pressure used in the N-S equation
        # to have the absolute pressure, one needs to have var("absolute_pressure")
    #    droplet_eqs+=IntegralObservableOutput("bottom_point")@"droplet_interface/droplet_axis"
        
        # Add all equations to the mesh domains
        self+=droplet_eqs@"droplet"+air_eqs@"air"
        self+=glass_eqs@"substrate"


    def presolve_temperature(self):
        self.default_mass_flux_coefficient=1*kilogram/(meter**2*second)
        # Start with a reasonable temperature profile by solving it before starting the simulation
        if not self.is_initialised():
            self.initialise()
        with self.select_dofs() as dofs:
            dofs.select("air/temperature","substrate/temperature","droplet/droplet_interface/masstrans_water") #,"droplet/droplet_interface/masstrans_water"
            dofs.select("droplet/droplet_interface/_lagr_conn_temperature_temperature","air/substrate_top/_lagr_conn_temperature_temperature")
            self.solve(max_newton_iterations=20)
        self.reapply_boundary_conditions()
        self.get_mesh("droplet").ignore_initial_condition=True
        self.get_mesh("air").ignore_initial_condition=True
        self.get_mesh("substrate").ignore_initial_condition=True
