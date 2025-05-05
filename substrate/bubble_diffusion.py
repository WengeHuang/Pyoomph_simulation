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
from pyoomph.equations.stokes_stream_func import *


# Create a plotter object
class Plotter(MatplotlibPlotter):
    def __init__(self, problem: Problem, filetrunk: str = "plot_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

    def define_plot(self):
        pr=cast("BubbleFreezingProblem",self.get_problem())
        h1=pr.get_ode("globals").get_value("ice_front_pos_y")
        h2=pr.channel_height
        r=pr.channel_width

        # Get pressure
        pres_expr=pr.get_ode("globals").get_value("gas_pressure")
        abs_pres=(pr.pressure_0-pr.absolute_pressure)
        pres_ratio=pres_expr/abs_pres
        pres_expanded=pr.get_ode("globals").get_code_gen().expand_placeholders(pres_ratio,True)
        pres=round(float(pres_expanded),3) # To float
        
        # Get mass
        mass=round(float(pr.get_ode("globals").get_value("gas_mass")/pr.initial_volume/pr.rho_air_init),3)
        
        # Get volume
        vol=round(float(pr.get_ode("globals").get_value("gas_volume")/pr.initial_volume),5)

        # Set view
        self.set_view(xmin=-r,xmax=r,ymin=h1-0.025*h2,ymax=1.2*(h1+h2))

        # Colorbars
        cbT=self.add_colorbar("c [kg/m$^3$]", position="top right", cmap="coolwarm")
        cbv=self.add_colorbar("velocity [$\mu$m/s]", position="top left", factor=1e6)

        # Add plots
        self.add_plot("liquid/c", colorbar=cbT)
        self.add_plot("liquid/velocity", colorbar=cbv, transform="mirror_x")
        self.add_plot("liquid/velocity", mode="arrows", transform="mirror_x")
        self.add_plot("liquid/ice_front",transform=[None,"mirror_x"],linecolor="purple")
        self.add_plot("liquid/bubble",transform=[None,"mirror_x"],linecolor="green")
        
        # Add text
        t=round(float(pr.get_current_time()/second),3)
        txt = self.add_text("t = {}\np = {}\nV = {}\nm = {}".format(t,pres,vol,mass), "top center", bbox=dict(facecolor='wheat', alpha=1, boxstyle="round"))
        txt.ymargin-=0.025


# Create a mesh
class Mesh(GmshTemplate):
    def __init__(self, bubble_height: float = 125*micro*meter, bubble_radius: float = 10*micro*meter, channel_height: float = 300*micro*meter, channel_width: float = 250*micro*meter):
        super().__init__()
        self.channel_height=channel_height
        self.bubble_height=bubble_height
        self.bubble_radius=bubble_radius
        self.channel_width=channel_width

    def define_geometry(self):
        
        # Create points
        p00, pW0= self.point(0,0, size=0.025*self.default_resolution), self.point(self.channel_width,0)
        p0H2, pWH2=self.point(0,self.channel_height,size=0.025*self.default_resolution), self.point(self.channel_width,self.channel_height,size=2*self.default_resolution)
        p0BH, pRBHB,p0BHplus,p0BHminus=self.point(0,self.bubble_height),self.point(self.bubble_radius,self.bubble_height, size=0.025*self.default_resolution),self.point(0,self.bubble_height+self.bubble_radius, size=0.025*self.default_resolution),self.point(0,self.bubble_height-self.bubble_radius, size=0.025*self.default_resolution)

        # Create lines
        self.create_lines(p0BHminus,"liquid_axis",p00,"ice_front",pW0,"liquid_wall",pWH2,"liquid_top",p0H2,"liquid_axis",p0BHplus)
        self.circle_arc(pRBHB,p0BHplus,center=p0BH,name="bubble")
        self.circle_arc(pRBHB,p0BHminus,center=p0BH,name="bubble")

        # Create surfaces
        self.plane_surface("liquid_axis","liquid_top","liquid_wall","ice_front","bubble",name="liquid")


# Create a problem
class BubbleFreezingProblem(Problem):
    def __init__(self):
        super().__init__()

        # Define the model, 5 models are going to be solved:
        # 1. Pure diffusion of uniform the air concentration (i.e. c0 concentration on top and bottom) in the liquid generates bubble growth, no moving ice front
        # 2. Pure diffusion of non-uniform air concentration in the liquid generates bubble growth, no moving ice front
        # 3. Pure diffusion of non-uniform air concentration in the liquid generates bubble growth, moving ice front
        # 4. Advection of concentration due to temperature in the liquid generates bubble growth, moving ice front
        # 5. Advection of concentration due to temperature (reversed MaT effects) in the liquid generates bubble growth, no moving ice front
        self.pure_diffusion_static_uniform_c_model = False
        self.pure_diffusion_static_model = False
        self.pure_diffusion_moving_model = False
        self.advection_temperature_model = False
        self.advection_fake_temperature_model = False

        # Define dimensions
        self.channel_height=400*micro*meter # height of the liquid
        self.channel_width=250*micro*meter # radius of the cylinders
        self.bubble_height=125*micro*meter # height of the bubble from ice front
        self.bubble_radius=10*micro*meter # radius of the bubble
        self.resolution=2.5 # mesh resolution

        # Define physical parameters of Ice and Liquid
        self.gravity=9.81*meter/second**2 * vector(0,-1) # gravity direction and strength       
        self.ice_front_speed=1e-5*meter/second # speed of the ice front
        self.rho_liq=(1000-0.232*var("c")/scale_factor("c"))*kilogram/meter**3  # density of liquid (temperature dependency is added later)
        self.mu_liq=1*milli*pascal*second # dynamic viscosity of liquid
        self.T0=0*celsius # temperature of ice front
        self.Tair=20*celsius # temperature of air
        self.bubble_T = self.T0 + self.bubble_height/self.channel_height*(self.Tair-self.T0) # temperature of the droplet
        self.K=0.03 # coefficient
        self.c0=0.01*kilogram/meter**3 # top boundary concentration
        self.csat=0.0286*kilogram/meter**3 # saturation concentration
        self.thermal_diffusivity = 1.4558e-7*meter**2/second

        # Bubble properties
        self.surface_tension=0.0721*newton/meter # surface tension (temperature dependency is added later)
        self.sigma0=0.0721*newton/meter # surface tension at c=0
        
        # Gas in bubble properties
        self.initial_volume=4/3*pi*self.bubble_radius**3 # initial volume of the bubble
        self.R_air=287.058*joule/(kilogram*kelvin) # gas constant of air
        self.rho_air_init=1.2754*kilogram/meter**3 # density of air
        self.diffusivity=1e-9*meter**2/second # diffusivity of air in water
        self.pressure_0=self.rho_air_init*self.R_air*self.bubble_T # pressure initially
        self.absolute_pressure=self.pressure_0-2*self.sigma0/self.bubble_radius # pressure initially

        # Values at current time step
        self.rho_air=var("gas_mass", domain="globals")/var("gas_volume", domain="globals") # density of air
        self.mass_transfer_rate = self.diffusivity*dot(var("normal"),grad(var("c", domain="liquid"))) # mass transfer rate

        # Molar mass of air
        self.soret=0.03/kelvin

        # Add plotter
        # self.plotter=[Plotter(self),ZoomedPlotter(self)]
        self.plotter=Plotter(self)

    # Presolve the gas phase
    def presolve_gas_phase(self,*,globally_convergent_newton:bool=False,max_newton_iterations:Union[None,int]=None):
        if not self.is_initialised():
            self.initialise()
        doflist:Set[str]=set()
        doflist.add("liquid/c")
        if self.advection_temperature_model or self.advection_fake_temperature_model:
            doflist.add("liquid/temperature")
        if len(doflist)==0:
            return
        with self.select_dofs() as dofs:
            dofs.select(*doflist)
            self.solve(do_not_set_IC=True,globally_convergent_newton=globally_convergent_newton,max_newton_iterations=max_newton_iterations)
            self.timestepper.set_num_unsteady_steps_done(0)
            self.get_mesh("liquid").ignore_initial_condition = True

    # Stop simulations when the bubble is too close to the ice front
    def actions_after_newton_solve(self):
        super().actions_after_newton_solve()
        y_pos = self.get_ode("globals").get_value("bubble_y")
        radius = self.get_ode("globals").get_value("bubble_radius")
        ice_front_y_pos = self.get_ode("globals").get_value("ice_front_pos_y")
        gas_volume = self.get_ode("globals").get_value("gas_volume")
        if float((y_pos-radius-ice_front_y_pos)/radius) < 0.05 or float((ice_front_y_pos+self.channel_height-radius-y_pos)/radius) < 0.05 or float(gas_volume/self.initial_volume) > 300:
            self.abort_current_run()

    # Define the problem
    def define_problem(self):
        
        # Redefine physical parameters according to the model
        self.absolute_pressure+=2*self.sigma0/self.bubble_radius # pressure initially
        if self.pure_diffusion_static_uniform_c_model:
            self.c0=0.04*kilogram/meter**3 # top boundary concentration
        if self.advection_temperature_model:
            TKelvin=var("temperature")/scale_factor("temperature")
            self.surface_tension+=-0.07275*0.002*(TKelvin-291.0)*newton/meter # surface tension
            self.sigma0+=-0.07275*0.002*(self.bubble_T/scale_factor("temperature")-291.0)*newton/meter # surface tension
            self.rho_liq+=((1000 * (TKelvin + 15.7914) / (508929 * (TKelvin - 205.02037)) * (TKelvin - 277.1363) ** 2))*kilogram/meter**3 # density of liquid
            chi = var("c")/self.rho_air_init
            self.add_equations(WeakContribution(self.rho_air_init*self.soret*self.diffusivity*grad(var("temperature"))*chi*(1-chi),grad(testfunction("c")))@"liquid")
        elif self.advection_fake_temperature_model:
            TKelvin=var("temperature")/scale_factor("temperature")
            self.surface_tension+=0.07275*0.002*(TKelvin-291.0)*newton/meter # surface tension
            self.sigma0+=0.07275*0.002*(self.bubble_T/scale_factor("temperature")-291.0)*newton/meter # surface tension
            self.rho_liq+=(1000 * ((TKelvin + 15.7914) / (508929 * (TKelvin - 205.02037)) * (TKelvin - 277.1363) ** 2))*kilogram/meter**3 # density of liquid
            chi = var("c")/self.rho_air_init
            self.add_equations(WeakContribution(self.rho_air_init*self.soret*self.diffusivity*grad(var("temperature"))*chi*(1-chi),grad(testfunction("c")))@"liquid")
        self.absolute_pressure+=-2*self.sigma0/self.bubble_radius # pressure initially
        self.curr_pres=var("gas_pressure", domain="globals")+self.absolute_pressure # current pressure

        # Set coordinates system
        self.set_coordinate_system("axisymmetric")

        # Scales for nondimensionalization
        self.set_scaling(spatial=self.bubble_radius,temporal=self.bubble_radius/self.ice_front_speed)
        self.set_scaling(velocity=self.ice_front_speed)
        self.set_scaling(pressure=self.pressure_0)
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
        if self.pure_diffusion_static_uniform_c_model or self.pure_diffusion_static_model:
            leqs+=DirichletBC(mesh_y=self.channel_height)@"liquid_top"
            leqs+=DirichletBC(mesh_y=0)@"ice_front"
        else:
            leqs+=IceFrontSpeed(ice_front_speed=-self.ice_front_speed)@"ice_front" # Advance the ice
            leqs+=IceFrontSpeed(ice_front_speed=self.ice_front_speed)@"liquid_top" # Mimic top far field as infinite
        leqs+=DirichletBC(mesh_x=self.channel_width)@"liquid_wall"

        # Describe advection diffusion equations for concentration
        leqs+=AdvectionDiffusionEquations(fieldnames="c",diffusivity=self.diffusivity, space="C2", wind=var("velocity") if (self.advection_temperature_model or self.advection_fake_temperature_model) else 0)
        leqs+=DirichletBC(c=self.c0)@"liquid_top"
        leqs+=DirichletBC(c=self.csat)@"bubble"
        if self.pure_diffusion_static_uniform_c_model:
            leqs+=DirichletBC(c=self.c0)@"ice_front"
        else:
            leqs+=DirichletBC(c=self.c0*(1+(1-self.K)/self.K))@"ice_front"
            leqs+=InitialCondition(c=self.c0*(1+(1-self.K)/self.K*exp(-self.ice_front_speed/self.diffusivity*(var("coordinate_y")))))

        # Add temperature effects
        if self.advection_temperature_model or self.advection_fake_temperature_model:
            leqs+=AdvectionDiffusionEquations("temperature", diffusivity=self.thermal_diffusivity, wind=var("velocity"))
            leqs+=DirichletBC(temperature=self.T0)@"ice_front"
            leqs+=DirichletBC(temperature=self.Tair)@"liquid_top"
            leqs+=StreamFunctionFromVelocity()+StreamFunctionFromVelocityInterface()@"bubble"+DirichletBC(streamfunc=0)@["liquid_top","ice_front","liquid_wall"]    

        # Mass transfer into the bubble and Marangoni effect
        leqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension, mass_transfer_rate=-self.mass_transfer_rate*self.rho_liq/self.rho_air, additional_normal_traction=var("gas_pressure", domain="globals"))@"bubble"

        # Gas volume and pressure
        leqs+=WeakContribution(self.mass_transfer_rate, testfunction("gas_mass", domain="globals"),dimensional_dx=True)@"bubble" # Mass transfer rate
        leqs+=WeakContribution(1/3*dot(var("mesh"),-var("normal")),testfunction("gas_volume", domain="globals"),dimensional_dx=True)@"bubble" # Volume of the bubble
        
        # Fix y position of bubble 
        leqs+=WeakContribution(var("mesh_y")-self.bubble_height,testfunction("fix_y_pos",domain="globals"))@"bubble"
        leqs+=WeakContribution(var("fix_y_pos",domain="globals"),testfunction("velocity_y"))

        #Lagrange multipliers
        geqs=ODEFileOutput()
        geqs+=GlobalLagrangeMultiplier(gas_mass=partial_t(var("gas_mass")))+Scaling(gas_mass=scale_factor("spatial")**3*self.rho_air_init)+TestScaling(gas_mass=scale_factor("temporal")/(self.rho_air_init*scale_factor("spatial")**3))
        geqs+=GlobalLagrangeMultiplier(gas_pressure=self.curr_pres)+Scaling(gas_pressure=scale_factor("pressure"))+TestScaling(gas_pressure=1/scale_factor("pressure"))
        geqs+=GlobalLagrangeMultiplier(gas_volume=-var("gas_volume"))+Scaling(gas_volume=scale_factor("spatial")**3)+TestScaling(gas_volume=1/scale_factor("spatial")**3)
        geqs+=GlobalLagrangeMultiplier(fix_y_pos=0)+Scaling(fix_y_pos=second**(-2)*meter**(-2)*kilogram)+TestScaling(fix_y_pos=1/scale_factor("spatial"))
        geqs+=WeakContribution(-self.rho_air*self.R_air*self.bubble_T, testfunction("gas_pressure")) # Ideal gas law
        geqs+=InitialCondition(gas_mass=self.initial_volume*self.rho_air_init)+InitialCondition(gas_volume=self.initial_volume)#+InitialCondition(gas_pressure=self.pressure_0-self.absolute_pressure)
        
        # Lagrange multipliers (not relevant for problem solving)
        geqs+=GlobalLagrangeMultiplier(ice_front_pos_y=0)+Scaling(ice_front_pos_y=scale_factor("spatial"))+TestScaling(ice_front_pos_y=1/scale_factor("spatial"))+InitialCondition(ice_front_pos_y=0)
        geqs+=GlobalLagrangeMultiplier(bubble_y=0)+Scaling(bubble_y=scale_factor("spatial"))+TestScaling(bubble_y=1/scale_factor("spatial"))+InitialCondition(bubble_y=self.bubble_height)
        geqs+=GlobalLagrangeMultiplier(bubble_radius=0)+Scaling(bubble_radius=scale_factor("spatial"))+TestScaling(bubble_radius=1/scale_factor("spatial"))+InitialCondition(bubble_radius=self.bubble_radius)
        geqs+=WeakContribution(var("bubble_radius")-square_root(var("gas_volume")/(4/3*pi),3),testfunction("bubble_radius"))
        leqs+=WeakContribution(var("ice_front_pos_y",domain="globals")-var("mesh_y"),testfunction("ice_front_pos_y",domain="globals"))@"ice_front"      
        leqs+=WeakContribution(var("bubble_y",domain="globals")-var("coordinate_y"),testfunction("bubble_y",domain="globals"))@"bubble"

        # Plotting
        leqs+=IntegralObservables(mass_transfer_rate=self.mass_transfer_rate)@"bubble"+IntegralObservables(volume_transfer_rate=self.mass_transfer_rate/self.rho_air)@"bubble"+IntegralObservables(volume=-1/3*dot(var("coordinate"),var("normal")))@"bubble"+IntegralObservables(volume_factor=-1/3*dot(var("coordinate"),var("normal"))/self.initial_volume)@"bubble"
        leqs+=IntegralObservableOutput("mass_transfer")@"bubble"
        leqs+=LocalExpressions(transfer_rate=self.mass_transfer_rate)@"bubble"
        leqs+=LocalExpressions(tangential_velocity=dot(var("velocity"),vector([-var("normal")[1],var("normal")[0]])))@"bubble"
        
        # Add equations
        self.add_equations(leqs@"liquid"+geqs@"globals")

with BubbleFreezingProblem() as p:
    p.plotter = Plotter(p)
    #p.plotter=None
    p.resolution=2.5
    p.pure_diffusion_static_uniform_c_model = False
    p.pure_diffusion_static_model = True
    p.pure_diffusion_moving_model = False
    p.advection_temperature_model = False
    p.advection_fake_temperature_model = False
    p.presolve_gas_phase()
    p.run(endtime=10*second, startstep=0.01*second, outstep=0.01*second, temporal_error=1)
    p.output_at_increased_time()