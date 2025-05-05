from pyoomph import *
from pyoomph.equations.multi_component import *
from pyoomph.equations.ALE import *
from pyoomph.utils.dropgeom import DropletGeometry
from pyoomph.equations.navier_stokes import *
from pyoomph.materials import *
from pyoomph.materials.mass_transfer import *
from pyoomph.meshes.remesher import Remesher2d
import pyoomph.materials.default_materials 
from pyoomph.meshes.zeta import *

from pyoomph.output.plotting import MatplotlibPlotter
from pyoomph.equations.contact_angle import *

from pyoomph.generic.problem import Problem

class SingleBubblePlotter(MatplotlibPlotter):
    def define_plot(self):
        pr=cast(SingleBubbleLaserProblem,self.get_problem())
        xrange=2*pr.R0
        self.background_color="darkgrey"
        self.set_view(-xrange,-0.3*xrange,xrange,xrange)
        cb_T=self.add_colorbar("temperature [°C]",offset=-273.15,position="bottom left")
        cb_u=self.add_colorbar("velocity [m/s]",position="bottom right",cmap="viridis")
        self.add_plot("liquid/temperature",colorbar=cb_T,transform="mirror_x")
        self.add_plot("bubble/temperature",colorbar=cb_T,transform="mirror_x")
        self.add_plot("substrate/temperature",colorbar=cb_T,transform="mirror_x")
        
        self.add_plot("liquid/velocity",colorbar=cb_u)
        self.add_plot("bubble/velocity",colorbar=cb_u)
        
        self.add_plot("bubble/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("liquid/velocity",mode="arrows",transform=[None,"mirror_x"])
        
        self.add_plot("liquid/interface",transform=[None,"mirror_x"])
        self.add_plot("liquid/liquid_substrate",transform=[None,"mirror_x"])
        self.add_plot("bubble/bubble_substrate",transform=[None,"mirror_x"])
        
        arrkey=self.add_arrow_key("bottom right",title="water mass transfer [kg/m²/s]")
        arrkey.ymargin=0.175
        arrkey.xmargin=0.15
        self.add_plot("liquid/interface/masstrans_water",arrowkey=arrkey,transform=[None,"mirror_x"])
        
        sb=self.add_scale_bar("bottom center")
        tl=self.add_time_label("bottom center")
        tl.unit="us"
        tl.yshift+=0.05
        

class SingleBubbleAxisymmMesh(GmshTemplate):
    def define_geometry(self):
        pr=cast(SingleBubbleLaserProblem,self.get_problem())
        geom=DropletGeometry(curv_radius=pr.R0,contact_angle=pr.theta0)
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
        
        p0S=self.point(0,-pr.Lhost,size=far_resolution)
        self.circle_arc(p0S,pL0,center=p00,name="substrate_far")
        self.line(p00,p0S,name="substrate_axis")
        self.plane_surface("substrate_axis","substrate_far","liquid_substrate","bubble_substrate",name="substrate")
        
        self.remesher=Remesher2d(self)
    
class SingleBubbleLaserProblem(Problem):
    def __init__(self):
        super().__init__()
        self.R0=4*micro*meter
        self.theta0=120*degree 
        self.typical_frequency=1*mega*hertz
        
        self.Lhost=100*micro*meter
        
        self.laser_power=12*milli*watt
        self.laser_conversion_efficiency=0.4
        self.laser_FWHM_width=10*micro*meter
        
        self.Troom=20*celsius
        self.Tbubble0=120*celsius
        
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

        self.contact_angle=120*degree # equilibrium contact angle        
        
    def get_laser_profile(self,r=var("coordinate_x")):
        w=self.laser_FWHM_width/square_root(2*log(2))
        E=exp(-r**2/w**2)
        return E**2 
        #return 0.1
        
    def define_problem(self):
        self.set_coordinate_system("axisymmetric")
        self.set_scaling(spatial=self.R0,temporal=1*micro*second,velocity=1*meter/second)
        self.liquid.set_reference_scaling_to_problem(self,temperature=100*celsius)

        
        if self.interface is None:
            self.interface=self.liquid | self.gas
            mass_transfer_model=HertzKnudsenSchrageMassTransferModel(self.liquid,self.gas)
            mass_transfer_model.projection_space="C2"
            #mass_transfer_model.sticking_coefficient=0
            self.interface.set_mass_transfer_model(mass_transfer_model)
            
        
        sigma0=self.liquid.evaluate_at_condition(self.interface.surface_tension,temperature=100*celsius)
        self.set_scaling(pressure=2*sigma0/self.R0)
        


            
        #self.define_named_var(absolute_pressure=1*atm+var("pressure"))
        self.define_named_var(absolute_pressure=1*atm + var("pressure"))
        self+=SingleBubbleAxisymmMesh()
        geom=DropletGeometry(curv_radius=self.R0,contact_angle=self.theta0)
        self.liquid.dynamic_viscosity=self.liquid.evaluate_at_condition("dynamic_viscosity",temperature=self.Troom)
        
        
        
        bubble_eqs=MeshFileOutput()
        liquid_eqs=MeshFileOutput()
        glass_eqs=MeshFileOutput()
        
        bubble_eqs+=CompositionFlowEquations(self.gas,isothermal=False,initial_temperature=self.Tbubble0,compo_space="C2",gravity=self.g*vector(0,-1))
        liquid_eqs+=CompositionFlowEquations(self.liquid,isothermal=False,initial_temperature=self.Troom,with_IC=False,compo_space="C2",gravity=self.g*vector(0,-1))
        glass_eqs+=TemperatureConductionEquation(self.substrate,space="C2")+InitialCondition(temperature=self.Troom)
        
        # Additional compressibility term for the temperature equation
        bubble_eqs+=WeakContribution(-(partial_t(var("pressure"))+dot(var("velocity"),grad(var("pressure")))),"temperature")
        
        R=square_root(dot(var("coordinate"),var("coordinate")))
        dist_factor=self.R0/R
        initial_T=self.Troom+minimum(1,dist_factor)*(self.Tbubble0-self.Troom)
        liquid_eqs+=InitialCondition(temperature=initial_T)
        glass_eqs+=InitialCondition(temperature=initial_T)
        
        bubble_eqs+=AxisymmetryBC()@"bubble_axis"
        liquid_eqs+=AxisymmetryBC()@"liquid_axis"
        glass_eqs+=AxisymmetryBC()@"substrate_axis"
        
        liquid_eqs+=LaplaceSmoothedMesh()
        bubble_eqs+=LaplaceSmoothedMesh()
        
        liquid_eqs+=DirichletBC(velocity_x=0)@"liquid_side"
        

        interf_eqs=MultiComponentNavierStokesInterface(self.interface)
        #interf_eqs=MultiComponentNavierStokesInterface(self.interface.surface_tension)
        interf_eqs+=ConnectFieldsAtInterface("temperature")+DirichletBC(_lagr_conn_temperature_temperature=0)@"liquid_substrate"
        #interf_eqs+=ConnectFieldsAtInterface("velocity")
        interf_eqs+=ConnectMeshAtInterface()



        
        liquid_eqs+=interf_eqs@"interface"
        #liquid_eqs+=MultiComponentNavierStokesInterface(self.interface.surface_tension)@"interface"
        
        bubble_eqs+=NeumannBC(temperature=-self.laser_power*self.laser_conversion_efficiency*self.get_laser_profile()/self.laser_FWHM_width**2)@"bubble_substrate"
        liquid_eqs+=NeumannBC(temperature=-self.laser_power*self.laser_conversion_efficiency*self.get_laser_profile()/self.laser_FWHM_width**2)@"liquid_substrate"
        
        glass_eqs+=TemperatureInfinityEquations(self.Troom)@"substrate_far"
        liquid_eqs+=TemperatureInfinityEquations(self.Troom)@"liquid_side"
        liquid_eqs+=TemperatureInfinityEquations(self.Troom)@"liquid_top"
        
        bubble_eqs+=ConnectFieldsAtInterface("temperature")@"bubble_substrate"
        bubble_eqs+=DirichletBC(_lagr_conn_temperature_temperature=0)@"bubble_substrate/interface"
        liquid_eqs+=ConnectFieldsAtInterface("temperature")@"liquid_substrate"
        
#contact angle
        #liquid_eqs+=NavierStokesContactAngle(self.contact_angle,wall_normal=vector(0,1),wall_tangent=vector(1,0))@"interface/liquid_substrate"
        #liquid_eqs+=NavierStokesContactAngle(self.contact_angle,wall_normal=vector(0,1))@"interface/liquid_substrate"
        
# contact angle must be defined on the liquid domain
# the wall_tangent vector (1, 0) and (-1, 0) makes the difference between ca and pi - ca        


        bubble_eqs+=NavierStokesSlipLength(self.sliplength)@"bubble_substrate"
        bubble_eqs+=DirichletBC(velocity_y=0,mesh_y=0)@"bubble_substrate"
        liquid_eqs+=NavierStokesSlipLength(self.sliplength)@"liquid_substrate"
        liquid_eqs+=DirichletBC(velocity_y=0,mesh_y=0)@"liquid_substrate"
        
        liquid_eqs+=RemeshWhen(self.remesh_options)
        bubble_eqs+=RemeshWhen(self.remesh_options)
            
        if self.pinned_contact_line:    # if pinned   
            # Pinned contact line means mesh_x is fixed.   
            # using mesh_x = True to fix it  
            liquid_eqs+=DirichletBC(mesh_x=True)@"liquid_substrate"
            bubble_eqs+=DirichletBC(mesh_x=True)@"bubble_substrate"

            from pyoomph.meshes.bcs import PythonDirichletBC
            unpin_contact_line_position=PythonDirichletBC(mesh_x=True)
            unpin_contact_line_position.unpin_instead=True
            liquid_eqs+=unpin_contact_line_position@"interface/liquid_substrate"
            unpin_contact_line_position=PythonDirichletBC(mesh_x=True)
            unpin_contact_line_position.unpin_instead=True
            bubble_eqs+=unpin_contact_line_position@"interface/liquid_substrate"
            liquid_eqs+=EnforcedBC(velocity_x=(var("mesh_x")-geom.base_radius)/scale_factor("temporal"))@"interface/liquid_substrate"
        else:            
            # Also let the substrate move
            glass_eqs+=LaplaceSmoothedMesh()+ElementSpace("C2")
            # Fix the height of the substrate (only move in x direction)
            glass_eqs+=DirichletBC(mesh_y=0)@"liquid_substrate"
            glass_eqs+=DirichletBC(mesh_y=0)@"bubble_substrate"
            # Fix the far field
            glass_eqs+=PinMeshCoordinates()@"substrate_far"
            # Connect the mesh at the liquid boundary
            
            liquid_eqs+=ConnectMeshAtInterface()@"liquid_substrate"
            # Connect the mesh at the gas boundary (but remove the connecting Lagrange multiplier at the corner to remove overconstraining)
            bubble_eqs+=(ConnectMeshAtInterface()+DirichletBC(_lagr_conn_mesh_x=True)@"interface")@"bubble_substrate" 
            
            #liquid_eqs+=DynamicContactLineEquations(UnpinnedContactLine(self.contact_angle,cl_speed_scale=None))@"interface/liquid_substrate"
            liquid_eqs+=NeumannBC(velocity_x=-self.interface.surface_tension*cos(self.contact_angle))@"interface/liquid_substrate"
            #liquid_eqs+=EnforcedBC(velocity_x=(var("normal_x")-cos(self.contact_angle))*scale_factor("velocity"))@"interface/liquid_substrate"
            liquid_eqs+=IntegralObservables(cartesian, contact_angle=atan2(-var("normal_y"),var("normal_x"))/degree)@"interface/liquid_substrate"
            liquid_eqs+=IntegralObservableOutput("ContactAngle")@"interface/liquid_substrate"
        
        
        bubble_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"interface"
        liquid_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"interface"
        glass_eqs+=AssignZetaCoordinatesByArclength(sort_along_axis="x+")@"substrate_far"
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
        
        bubble_eqs+=IntegralObservables(volume=1)
        bubble_eqs+=IntegralObservableOutput("bubble_evolution")
        #bubble_eqs+=IntegralObservables(volume=1,mass=self.gas.mass_density)
        self+=bubble_eqs@"bubble"+liquid_eqs@"liquid"
        self+=glass_eqs@"substrate"
        #self.add_equations
    


    def presolve_temperature(self):
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
        self.output()



problem=SingleBubbleLaserProblem()
problem+=SingleBubblePlotter(problem)
#problem.liquid.sample_all_properties_to_text_files(problem.get_output_directory("sample_liquid"),temperature=100*celsius)
T=1/problem.typical_frequency
dt=0.1*T
#problem.gas.mass_density=1*kilogram/meter**3
#problem.gas.mass_density=problem.gas.evaluate_at_condition("mass_density",temperature=100*celsius,absolute_pressure=1*atm)


if False:
    problem.interface=problem.liquid | problem.gas
    problem.interface.surface_tension=72*milli*newton/meter
    mass_transfer_model=HertzKnudsenSchrageMassTransferModel(problem.liquid,problem.gas)
    mass_transfer_model.projection_space="C2"
    mass_transfer_model.sticking_coefficient=0
    problem.interface.set_mass_transfer_model(mass_transfer_model)

if False:
    lvl=4
    problem+=RefineToLevel(lvl)@"liquid/interface/liquid_substrate"
    problem+=RefineToLevel(lvl)@"bubble/interface/bubble_substrate"
    problem+=RefineToLevel(lvl)@"substrate/liquid_substrate/bubble_substrate"

problem.presolve_temperature()
#problem.run(15*T,outstep=True,startstep=0.001*dt,temporal_error=1)
problem.run(15*T,outstep=dt,temporal_error=1)
#problem.output_at_increased_time()

