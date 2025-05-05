import sys
sys.path.append('../../problem_def')
from pyoomph import *
from pyoomph.generic.problem import Problem
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.output.plotting import *
from pyoomph.typings import List, Optional, Union
from pyoomph.expressions.units import *
from pyoomph.equations.advection_diffusion import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.ALE import *
from pyoomph.meshes.remesher import *
from pyoomph.materials.default_materials import *
from pyoomph.equations.poisson import *
from problem_def.problem import *
#from plotter import *
#from equations import *
#from mesh import *

# Create a plotter object
class Plotter(MatplotlibPlotter):
    def __init__(self, problem: Problem, filetrunk: str = "plot_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

    def define_plot(self):
        pr=cast("BubbleFreezingProblem",self.get_problem())
        h1=pr.get_ode("globals").get_value("ice_front_pos_y")
        h2=pr.channel_height
        r=pr.channel_width
        

        # Set view
        self.set_view(xmin=-r,xmax=r,ymin=h1-0.025*h2,ymax=1.2*(h1+h2))

        # Colorbars
        cbT=self.add_colorbar("c [kg/m$^3$]", position="top right", cmap="jet")
        cbv=self.add_colorbar("velocity [$\\mathrm{{\\mu}}$m/s]", position="top left", cmap="coolwarm", factor=1e6)

        # Add plots
        self.add_plot("liquid/c", colorbar=cbT)
        self.add_plot("liquid/velocity", colorbar=cbv, transform="mirror_x")
        self.add_plot("liquid/velocity", mode="arrows", transform="mirror_x")
        self.add_plot("liquid/ice_front",transform=[None,"mirror_x"],linecolor="purple")
        self.add_plot("liquid/bubble",transform=[None,"mirror_x"],linecolor="green")
        
        # Add text
        t=round(float(pr.get_current_time()/second),3)
        txt = self.add_text("t = {}".format(t), "top center", bbox=dict(facecolor='wheat', alpha=1, boxstyle="round"))
        txt.ymargin-=0.025


# Create a plotter object
class ZoomedPlotter(MatplotlibPlotter):
    def __init__(self, problem: Problem, filetrunk: str = "zoomed_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

    def define_plot(self):
        pr=cast("BubbleFreezingProblem",self.get_problem())

        # Set view
        xmin=-1.5*10*micro*meter
        xmax=1.5*10*micro*meter
        ymin=-1.5*10*micro*meter+125*micro*meter
        ymax=1.5*10*micro*meter+125*micro*meter
        self.set_view(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)

        # Colorbars
        cbT=self.add_colorbar("c [kg/m$^3$]", position="top right", cmap="jet")
        cbv=self.add_colorbar("velocity [$\mu$m/s]", position="top left", factor=1e6)

        # Add plots
        self.add_plot("liquid")
        self.add_plot("liquid/c", colorbar=cbT)
        self.add_plot("liquid/velocity", colorbar=cbv, transform="mirror_x")
        self.add_plot("liquid/velocity", mode="arrows", transform="mirror_x")
        self.add_plot("liquid/ice_front",transform=[None,"mirror_x"],linecolor="purple")
        self.add_plot("liquid/bubble",transform=[None,"mirror_x"],linecolor="green")

        # Add time label
        self.add_time_label("top center")

# Create a problem
class BubbleFreezingProblem(Problem):
    def __init__(self):
        super().__init__()

        # Define dimensions
        self.channel_height=300*micro*meter # height of the liquid
        self.channel_width=250*micro*meter # radius of the cylinders
        self.bubble_height=125*micro*meter # height of the bubble from ice front
        self.bubble_radius=10*micro*meter # radius of the bubble
        self.resolution=2.5 # mesh resolution The smaller it is the fineer the mesh would be

        # Define physical parameters of Ice and Liquid
        self.gravity=9.81*meter/second**2 * vector(0,-1) # gravity direction and strength       
        self.ice_front_speed=1e-5*meter/second # speed of the ice front
        self.rho_liq=(1000-0.232*var("c")/scale_factor("c"))*kilogram/meter**3  # density of liquid (temperature dependency is added later)
        self.mu_liq=1*milli*pascal*second # dynamic viscosity of liquid
        self.T0=0*celsius # temperature of ice front
        self.Tair=50*celsius # temperature of air
        self.K=0.03 # coefficient
        self.c0=0.01*kilogram/meter**3 # top boundary concentration
        self.csat=0.0286*kilogram/meter**3 # saturation concentration
        self.thermal_diffusivity = 1.4558e-7*meter**2/second

        # Bubble properties
        self.surface_tension=0.0721*newton/meter # surface tension (temperature dependency is added later)
        self.sigma0=0.0721*newton/meter # surface tension at c=0
        
        # Gas in bubble properties
        self.initial_volume=4/3*pi*self.bubble_radius**3 # initial volume of the bubble
        self.diffusivity=1e-9*meter**2/second # diffusivity of air in water
        self.R_air=287.058*joule/(kilogram*kelvin) # gas constant of air
        self.rho_air_init=1.2754*kilogram/meter**3 # density of air
        self.mu_air=1.729e-5*pascal*second # dynamic viscosity of air
        self.diffusivity=1e-9*meter**2/second # diffusivity of air in water
        self.pressure_0=self.rho_air_init*self.R_air*self.T0 # pressure initially
        self.absolute_pressure=self.pressure_0-2*self.sigma0/self.bubble_radius # absolute pressure

    # Presolve the gas phase
    def presolve_gas_phase(self,*,globally_convergent_newton:bool=False,max_newton_iterations:Union[None,int]=None):
        if not self.is_initialised():
            self.initialise()
        doflist:Set[str]=set()
        doflist.add("liquid/c")
        doflist.add("liquid/temperature")
        if len(doflist)==0:
            return
        with self.select_dofs() as dofs:
            dofs.select(*doflist)
            self.solve(do_not_set_IC=True,globally_convergent_newton=globally_convergent_newton,max_newton_iterations=max_newton_iterations)
            self.timestepper.set_num_unsteady_steps_done(0)
            self.get_mesh("liquid").ignore_initial_condition = True

    # Define the problem
    def define_problem(self):
        
        # Redefine physical parameters according to the model
        TKelvin=var("temperature")/scale_factor("temperature")
        self.surface_tension+=0.07275*0.002*(TKelvin-291.0)*newton/meter # surface tension
        self.sigma0+=0.07275*0.002*(self.Tair/scale_factor("temperature")-291.0)*newton/meter # surface tension
        self.rho_liq+=(1000 * ((TKelvin + 15.7914) / (508929 * (TKelvin - 205.02037)) * (TKelvin - 277.1363) ** 2))*kilogram/meter**3 # density of liquid
    
        # Add plotter
        self.plotter=[Plotter(self),ZoomedPlotter(self)]

        # Set coordinates system
        self.set_coordinate_system("axisymmetric")

        # Scales for nondimensionalization
        self.set_scaling(spatial=self.bubble_radius,temporal=self.bubble_radius/self.ice_front_speed)
        self.set_scaling(velocity=self.ice_front_speed)
        self.set_scaling(pressure=self.sigma0/scale_factor("spatial"))
        self.set_scaling(c=self.c0)
        self.set_scaling(temperature=self.T0)
        
        # Add mesh
        mesh=Mesh(bubble_radius=self.bubble_radius, bubble_height=self.bubble_height, channel_height=self.channel_height, channel_width=self.channel_width)
        mesh.default_resolution=self.resolution
        mesh.remesher=Remesher2d(mesh)
        mesh.mesh_mode="tris"
        self.add_mesh(mesh)

         # Mesh and output stuff
        leqs=MeshFileOutput()
        leqs+=TextFileOutput()
        leqs+=TextFileOutput()@"bubble"
        leqs+=PseudoElasticMesh()
        leqs+=RemeshWhen(RemeshingOptions())

        # Describe velocity equations for liquid
        leqs+=NavierStokesEquations(mass_density=self.rho_liq,dynamic_viscosity=self.mu_liq,gravity=self.gravity, boussinesq=True, mode="TH")
        leqs+=NoSlipBC()@"liquid_wall"
        leqs+=AxisymmetryBC()@"liquid_axis"
        leqs+=NoSlipBC()@"ice_front" # No slip condition at the ice front
        leqs+=IceFrontSpeed(ice_front_speed=-self.ice_front_speed)@"ice_front" # Advance the ice
        leqs+=IceFrontSpeed(ice_front_speed=self.ice_front_speed)@"liquid_top" # Mimic top far field as infinite
        leqs+=DirichletBC(mesh_x=self.channel_width)@"liquid_wall"

        # Describe advection diffusion equations for concentration
        # leqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C1", wind=0)
        leqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C1", wind=var("velocity"))
        leqs+=AdvectionDiffusionInfinity(c=self.c0)@"liquid_top"
        leqs+=DirichletBC(c=self.csat)@"bubble"
        leqs+=DirichletBC(c=self.c0*(1+(1-self.K)/self.K))@"ice_front"
        leqs+=InitialCondition(c=self.c0*(1+(1-self.K)/self.K*exp(-self.ice_front_speed/self.diffusivity*(var("coordinate_y")))))

        # Add temperature effects
        leqs+=AdvectionDiffusionEquations("temperature", diffusivity=self.thermal_diffusivity, space="C1", wind=var("velocity"))
        leqs+=AdvectionDiffusionInfinity(temperature=self.Tair)@"liquid_top"
        leqs+=DirichletBC(temperature=self.T0)@"ice_front"

        # Static bubble model
        leqs+=DirichletBC(mesh_x=True, mesh_y=True)@"bubble"

        # Mass transfer into the bubble and Marangoni effect
        leqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension, static_interface=True, additional_normal_traction=-self.pressure_0)@"bubble"

        #Lagrange multipliers
        geqs=ODEFileOutput()
        geqs+=GlobalLagrangeMultiplier(ice_front_pos_y=0)+Scaling(ice_front_pos_y=scale_factor("spatial"))+TestScaling(ice_front_pos_y=1/scale_factor("spatial"))+InitialCondition(ice_front_pos_y=0)
        leqs+=WeakContribution(var("ice_front_pos_y",domain="globals")-var("mesh_y"),testfunction("ice_front_pos_y",domain="globals"))@"ice_front"      

        # Plotting
        leqs+=LocalExpressions(tangential_velocity=dot(var("velocity"),vector([-var("normal")[1],var("normal")[0]])))@"bubble"
        
        # Add equations
        self.add_equations(leqs@"liquid"+geqs@"globals")



with BubbleFreezingProblem() as p:
    p.presolve_gas_phase()
    p.run(10*second, startstep=0.01*second, maxstep=0.2*second, temporal_error=1)
    p.output_at_increased_time()