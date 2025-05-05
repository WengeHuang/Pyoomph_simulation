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


class SingleBubblePlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr = cast(SingleBubbleLaserProblem, self.get_problem())
        xrange = 3*pr.R0 # view range
        self.background_color = "darkgrey"
        self.set_view(-3*xrange, -0.5*xrange, 3*xrange, 3*xrange)  # -x_min, -ymin, x_max, y_max of the view window
        cb_T = self.add_colorbar("temperature [°C]", offset=-273.15, position = "bottom left")
        cb_u = self.add_colorbar("velocity [m/s]", position = "bottom right", cmap = "viridis")
        self.add_plot("liquid/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("bubble/temperature", colorbar=cb_T, transform = "mirror_x")
        self.add_plot("substrate/temperature",colorbar=cb_T,transform="mirror_x")
        
        self.add_plot("liquid/velocity",colorbar=cb_u)
        self.add_plot("bubble/velocity",colorbar=cb_u)
        
        self.add_plot("bubble/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("liquid/velocity",mode="arrows",transform=[None,"mirror_x"])
        
        self.add_plot("liquid/interface",transform=[None,"mirror_x"])
        self.add_plot("liquid/liquid_substrate",transform=[None,"mirror_x"])
        self.add_plot("bubble/bubble_substrate",transform=[None,"mirror_x"])

        arrkey=self.add_arrow_key("bottom right",title="water mass transfer [kg/m²/s]")

        arrkey.rangemode="fixed"
        arrkey.set_range(10)


        arrkey.ymargin=0.175
        arrkey.xmargin=0.15
        self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform=[None,"mirror_x"])
        
        sb=self.add_scale_bar("bottom center")
        #tl=self.add_time_label("bottom center")
        #tl.unit="us"
        #tl.yshift+=0.05

        #alph=self.add_text("alpha={:5g}".format(pr.alpha.value),position="bottom center",bbox=dict(boxstyle='round', facecolor='wheat', alpha=1),textsize=16)
        #alph.yshift+=0.075

# Mesh class. An axisymmetric bubble with liquid around on a substrate
class SingleBubbleAxisymmMesh(GmshTemplate):
    def define_geometry(self):
        pr=cast(SingleBubbleLaserProblem,self.get_problem())
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
        
        #p0S=self.point(0,-pr.Lhost,size=far_resolution)
        #self.circle_arc(p0S,pL0,center=p00,name="substrate_far")
        #self.line(p00,p0S,name="substrate_axis")
        #self.plane_surface("substrate_axis","substrate_far","liquid_substrate","bubble_substrate",name="substrate")
        p0S=self.point(0,-pr.thick)
        pLS=self.point(pr.Lhost,-pr.thick,size=far_resolution)
        self.create_lines(p00,"substrate_axis",p0S,"substrate_base",pLS,"substrate_side",pL0)
        self.plane_surface("substrate_axis","substrate_base","substrate_side","liquid_substrate","bubble_substrate",name="substrate")
        # attach a remesher for meshing
        self.remesher = Remesher2d(self)


# Problem class. Combining the equations, meshes, etc to a full problem
class SingleBubbleLaserProblem(Problem):
    def __init__(self):
        super().__init__()

        # Default parameters
        self.R0=4*micro*meter
        self.theta0=120*degree 
        self.typical_frequency=1*mega*hertz
        
        self.Lhost=100*micro*meter
        self.thick=10*micro*meter
        self.sticking_coefficient=1

        self.laser_power=12*milli*watt
        self.laser_conversion_efficiency=0.4
        self.laser_FWHM_width=10*micro*meter
        
        self.Troom=20*celsius
        self.Tbubble0=120*celsius
        self.Tsubstrate=120*celsius
        self.gas=get_pure_gas("water")
        self.liquid=get_pure_liquid("water")
        self.substrate=get_pure_solid("borosilicate")
        
        #self.pinned_contact_line=True
        self.pinned_contact_line=False
        
        self.sliplength=0.01*micro*meter
        
        self.interface=None
        
        self.g=9.81*meter/second**2
        
        self.remesh_options=RemeshingOptions()
        self.remesh_options.active=True  

        self.alpha = self.define_global_parameter(alpha = 0.1) 

    def get_interface_properties(self):
        if self.interface is None:
            self.interface=self.liquid | self.gas
            mass_transfer_model=HertzKnudsenSchrageMassTransferModel(self.liquid,self.gas)
            mass_transfer_model.projection_space="C2"
            mass_transfer_model.sticking_coefficient=self.alpha
            mass_transfer_model.with_prosperetti_term=True
            self.interface.set_mass_transfer_model(mass_transfer_model)
        return self.interface

    def get_laser_profile(self,r=var("coordinate_x")):
        w=self.laser_FWHM_width/square_root(2*log(2))
        E=exp(-r**2/w**2)
        return E**2 
        #return 0.1

    def define_problem(self):
        # Switch to an axisymmetric coordinate system
        self.set_coordinate_system("axisymmetric")
        # Select some reasonable scales
        self.set_scaling(spatial=self.R0,temporal=1*micro*second,velocity=1*meter/second)
        # And set the remaining scales by the liquid properties
        self.liquid.set_reference_scaling_to_problem(self,temperature=100*celsius)

        # The pressure scale is given by the Laplace pressure, so we must evaluate the surface tension at some constant temperature
        self.get_interface_properties()                    
        sigma0=self.liquid.evaluate_at_condition(self.interface.surface_tension,temperature=100*celsius)
        self.set_scaling(pressure=2*sigma0/self.R0)

        # the absolute pressure is the ambient pressure plus the pressure due to the Navier-Stokes
        # absolute_pressure is used e.g. in the ideal gas law for self.gas.mass_density
        self.define_named_var(absolute_pressure=1*atm + var("pressure"))

        # Add the mesh
        self+=SingleBubbleAxisymmMesh()
        geom=DropletGeometry(curv_radius=self.R0,contact_angle=self.contact_angle)     
        # Simplify: Use a constant viscosity (temperature and pressure independent)
        ## if we have something temperature and pressure dependent, how should we do
        self.liquid.dynamic_viscosity=self.liquid.evaluate_at_condition("dynamic_viscosity",temperature=self.Troom)  

        # Begin to assemble the equations for each of the three domains
        bubble_eqs=MeshFileOutput() # Output. PVD/VTU files can be loaded with Paraview
        liquid_eqs=MeshFileOutput()
        glass_eqs=MeshFileOutput()

        # Bubble is a non-isothermal, compressible Navier-Stokes equation (and potentially also with composition dynamics in case of gas mixtures)        
        ## how the compressibility is added
        bubble_eqs+=CompositionFlowEquations(self.gas,isothermal=False,initial_temperature=self.Tbubble0,compo_space="C2",gravity=self.g*vector(0,-1))       
        bubble_eqs+=InitialCondition(pressure=2*sigma0/self.R0) # Start with some reasonable pressure, since the density depends on it
        # Liquid is about the same, just Navier-Stokes
        liquid_eqs+=CompositionFlowEquations(self.liquid,isothermal=False,initial_temperature=self.Troom,with_IC=False,compo_space="C2",gravity=self.g*vector(0,-1))
        # Substrate just a conduction equation for the temperature. This is included in the CompositionFlowEquations as well, since we set isothermal=False there
        glass_eqs+=TemperatureConductionEquation(self.substrate,space="C2")+InitialCondition(temperature=self.Troom)

        # Additional compressibility term for the temperature equation
        bubble_eqs+=WeakContribution(-(partial_t(var("pressure"))+dot(var("velocity"),grad(var("pressure")))),"temperature")

        # Get some guess for the temperature profile
        R=square_root(dot(var("coordinate"),var("coordinate")))
        dist_factor=self.R0/R
        initial_T=self.Troom+minimum(1,dist_factor)*(self.Tbubble0-self.Troom)
        liquid_eqs+=InitialCondition(temperature=initial_T)
        glass_eqs+=InitialCondition(temperature=initial_T)

        # Simple moving mesh dynamics
        liquid_eqs+=LaplaceSmoothedMesh()
        bubble_eqs+=LaplaceSmoothedMesh()

         # Axis of symmetry, sets e.g. velocity_x=0
        bubble_eqs+=AxisymmetryBC()@"bubble_axis"
        liquid_eqs+=AxisymmetryBC()@"liquid_axis"
        glass_eqs+=AxisymmetryBC()@"substrate_axis"

        # Far field stuff. Must be set to 0 or the hydrostatic pressure must be balanced at the side
        liquid_eqs+=DirichletBC(velocity_x=0)@"liquid_side"   

        # Interface: Surface tension, mass transfer, velocity connection, Stefan flow, Marangoni flow, etc.
        interf_eqs=MultiComponentNavierStokesInterface(self.interface) 
         # Temperature must be connected (but the connection must be deactivated at the contact line, will be enforced via the substrate connection)       
        interf_eqs+=ConnectFieldsAtInterface("temperature")+DirichletBC(_lagr_conn_temperature_temperature=0)@"liquid_substrate"
        #interf_eqs+=ConnectMeshAtInterface()
        interf_eqs+=ConnectMeshAtInterface() + DirichletBC(_lagr_conn_mesh_x=0)@"liquid_substrate" # Mesh must stay connected        
        liquid_eqs+=interf_eqs@"interface"         

        # Heating by the laser
        #bubble_eqs+=NeumannBC(temperature=-self.laser_power*self.laser_conversion_efficiency*self.get_laser_profile()/self.laser_FWHM_width**2)@"bubble_substrate"
        #liquid_eqs+=NeumannBC(temperature=-self.laser_power*self.laser_conversion_efficiency*self.get_laser_profile()/self.laser_FWHM_width**2)@"liquid_substrate"
        glass_eqs+=DirichletBC(temperature=self.Tsubstrate)@"substrate_base"
        # Far field conditions for the temperature
        glass_eqs+=TemperatureInfinityEquations(self.Troom)@"substrate_side"
        liquid_eqs+=TemperatureInfinityEquations(self.Troom)@"liquid_side"
        liquid_eqs+=TemperatureInfinityEquations(self.Troom)@"liquid_top"

        # Connecting the temperature also to the substrate
        bubble_eqs+=ConnectFieldsAtInterface("temperature")@"bubble_substrate"
        # deactivated at the contact line
        bubble_eqs+=DirichletBC(_lagr_conn_temperature_temperature=0)@"bubble_substrate/interface"
        liquid_eqs+=ConnectFieldsAtInterface("temperature")@"liquid_substrate"

        # Slip length and further conditions
        bubble_eqs+=NavierStokesSlipLength(self.sliplength)@"bubble_substrate"
        bubble_eqs+=DirichletBC(velocity_y=0,mesh_y=0)@"bubble_substrate"
        liquid_eqs+=NavierStokesSlipLength(self.sliplength)@"liquid_substrate"
        liquid_eqs+=DirichletBC(velocity_y=0,mesh_y=0)@"liquid_substrate"

        if self.pinned_contact_line:
            # Pinned contact line: We can pin the mesh x position in  this case
            liquid_eqs+=DirichletBC(mesh_x=True)@"liquid_substrate"
            bubble_eqs+=DirichletBC(mesh_x=True)@"bubble_substrate"
            # But we must unpin it at the contact line directly
            from pyoomph.meshes.bcs import PythonDirichletBC
            unpin_contact_line_position=PythonDirichletBC(mesh_x=True)
            unpin_contact_line_position.unpin_instead=True
            liquid_eqs+=unpin_contact_line_position@"interface/liquid_substrate"
            unpin_contact_line_position=PythonDirichletBC(mesh_x=True)
            unpin_contact_line_position.unpin_instead=True
            bubble_eqs+=unpin_contact_line_position@"interface/liquid_substrate"
            liquid_eqs+=EnforcedBC(velocity_x=(var("mesh_x")-geom.base_radius)/scale_factor("temporal"))@"interface/liquid_substrate"
        else:            
            # Also let the substrate mesh move
            glass_eqs+=LaplaceSmoothedMesh()+ElementSpace("C2")
            # Fix the height of the substrate (only move in x direction)
            glass_eqs+=DirichletBC(mesh_y=0)@"liquid_substrate"
            glass_eqs+=DirichletBC(mesh_y=0)@"bubble_substrate"
            # Fix the far field
            glass_eqs+=PinMeshCoordinates()@"substrate_side"
            # Connect the mesh at the liquid boundary            
            liquid_eqs+=ConnectMeshAtInterface()@"liquid_substrate"
            # Connect the mesh at the gas boundary (but remove the connecting Lagrange multiplier at the corner to remove overconstraining)
            bubble_eqs+=(ConnectMeshAtInterface()+DirichletBC(_lagr_conn_mesh_x=True)@"interface")@"bubble_substrate" 
            # Neumann force to obtain the desired contact angle
            liquid_eqs+=NeumannBC(velocity_x=-self.interface.surface_tension*cos(self.contact_angle))@"interface/liquid_substrate"

        # Get the initial volume
        #V0=DropletGeometry(curv_radius=self.R0,contact_angle=self.theta0).volume
        # add a single Lagrange multiplier for the additional pressure required in the bubble to maintain the volume
        # We subtract -V0*pBadd_test from the equation of this Lagrange multiplier (equation_contribution=-V0)
        # pBadd should have and inverse scaling of the testfunction of the pressure, so that (pBadd,testfunction("pressure")) is nondimensional (scaling=...)
        # and we nondimensionalize pBadd_test=testfunction("pBadd") by 1/V0, so that the volume constraint is nondimensional (testscaling=1/V0)
        #pBadd,pBadd_test=self.add_global_dof("pBadd",equation_contribution=-V0,scaling=1/test_scale_factor("pressure"),testscaling=1/V0)
        # But we still must add the actual volume to the constraint. So we add (1,pBadd_test)_bubble to the residuals.
        # It must be integrated with dimensions, so that we get a volume contribution in m^3. Since pBadd_test scales with 1/V0, it is nondimensional
        # In total, pBadd will be adjusted so that the constraint (1,pBadd_test)_bubble-V0*pBadd_test=0 is satisfied.
        #bubble_eqs+=WeakContribution(1,pBadd_test,dimensional_dx=True)
        # Finally, we add this to the pressure of the bubble. This will adjust the pressure and lead to expansion/shrinkage of the bubble until V=V0
        #bubble_eqs+=WeakContribution(pBadd,"pressure")
        #bubble_eqs+=EnforceVolumeByPressure(V0)


        # Remesh when necessary
        liquid_eqs+=RemeshWhen(self.remesh_options)
        bubble_eqs+=RemeshWhen(self.remesh_options)
        # Zeta interpolation upon remeshing
        bubble_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"interface"
        liquid_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"interface"
        glass_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"substrate_side"
        liquid_eqs+=AssignZetaCoordinatesByEulerianCoordinate("x")@"liquid_substrate"
        bubble_eqs+=AssignZetaCoordinatesByEulerianCoordinate("x")@"bubble_substrate"
        liquid_eqs+=AssignZetaCoordinatesByEulerianCoordinate("y")@"liquid_axis"
        bubble_eqs+=AssignZetaCoordinatesByEulerianCoordinate("y")@"bubble_axis"
        glass_eqs+=AssignZetaCoordinatesByEulerianCoordinate("y")@"substrate_axis"
        glass_eqs+=AssignZetaCoordinatesByEulerianCoordinate("x")@"liquid_substrate"
        glass_eqs+=AssignZetaCoordinatesByEulerianCoordinate("x")@"bubble_substrate"
        bubble_eqs+=AssignZetaCoordinatesByEulerianCoordinate("y")@"bubble_axis"
        liquid_eqs+=AssignZetaCoordinatesByEulerianCoordinate("y")@"liquid_side"
        liquid_eqs+=AssignZetaCoordinatesByEulerianCoordinate("x")@"liquid_top"
        
        # Measure the volume (integral over 1 ) and the mass (integral over rho)
        bubble_eqs+=IntegralObservables(volume=1,mass=self.gas.mass_density)
        bubble_eqs+=IntegralObservableOutput("bubble_evolution") # write it to bubble_evolution.txt
        
        # Add all equations to the mesh domains
        self+=bubble_eqs@"bubble"+liquid_eqs@"liquid"
        self+=glass_eqs@"substrate"


    def presolve_temperature(self):
        # Start with a reasonable temperature profile by solving it before starting the simulation
        if not self.is_initialised():
            self.initialise()
        with self.select_dofs() as dofs:
            dofs.select("liquid/temperature","bubble/temperature","substrate/temperature","liquid/interface/masstrans_water")
            dofs.select("liquid/interface/_lagr_conn_temperature_temperature","bubble/bubble_substrate/_lagr_conn_temperature_temperature","liquid/liquid_substrate/_lagr_conn_temperature_temperature")
            self.solve()
        self.reapply_boundary_conditions()
        self.get_mesh("liquid").ignore_initial_condition=True
        self.get_mesh("bubble").ignore_initial_condition=True
        self.get_mesh("substrate").ignore_initial_condition=True