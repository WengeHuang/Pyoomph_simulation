# Import pyoomph and all used equations and utilities
from pyoomph import * # basic stuff
 
from pyoomph.equations.advection_diffusion import *
#from pyoomph.equations.multi_component import * #multi-component flow equations
from pyoomph.equations.ALE import * # moving mesh equations
from pyoomph.equations.navier_stokes import * #Navier-Stokes stuff, e.g. slip length
from pyoomph.equations.contact_angle import * #Contact angle equations
from pyoomph.equations.poisson import * #Navier-Stokes stuff, e.g. slip length

from pyoomph.utils.dropgeom import DropletGeometry # droplet geometry calculations
 
from pyoomph.materials import * # Material API
from pyoomph.materials.mass_transfer import * # Mass transfer models
import pyoomph.materials.default_materials # Default materials like water, water vapor, glass
 
from pyoomph.meshes.remesher import Remesher2d # Remeshing on larger deformation
from pyoomph.meshes.zeta import * # Zeta coordinate interpolation upon remeshing
 
from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools

 
# Plotter calss. Will generate figures
class SingleBubblePlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr = cast(SingleBubbleProblem, self.get_problem())
        xrange = 3*pr.R0 # view range
        self.background_color = "darkgrey"
        self.set_view(-xrange, -0.3*xrange, xrange, xrange)  # -x_min, -ymin, x_max, y_max of the view window
        #self.set_view(-0.5*xrange, -0.05*xrange, 0, 0.4*xrange)
        cb_T = self.add_colorbar("temperature [°C]",  position = "top right")
        cb_c = self.add_colorbar("velocity [m/s]", position = "bottom right", cmap = "viridis")
        #cb_c=self.add_colorbar("c [g/m^3]",cmap="Blues",position="top left")
        #showing mesh
        #self.add_plot("liquid")
        #self.add_plot("Pt")
        #self.add_plot("substrate")
 
 
        #self.add_plot("liquid/c", colorbar=cb_c, transform = "mirror_x")
        #self.add_plot("bubble/c", colorbar=cb_c, transform = "mirror_x")
        #self.add_plot("Pt/J", colorbar=cb_c,transform=[None,"mirror_x"])
        #self.add_plot("substrate/temperature", colorbar=cb_T)
        self.add_plot("liquid/J",colorbar=cb_T,transform="mirror_x")
       
        #self.add_plot("liquid/cc",colorbar=cb_T)
        #self.add_plot("bubble/velocity",colorbar=cb_u)
        #self.add_plot("bubble/c",colorbar=cb_u)
       
        #self.add_plot("bubble/velocity",mode="arrows",transform=[None,"mirror_x"])
        self.add_plot("liquid/J",mode="streamlines",transform="mirror_x")
        #self.add_plot("Pt/J",mode="streamlines",transform="mirror_x")
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
        tl=self.add_time_label("bottom center")
 
# Mesh class. An axisymmetric bubble with liquid around on a substrate
class SingleBubbleAxisymmMesh_angle(GmshTemplate):
    def define_geometry(self):
        pr=cast(SingleBubbleProblem,self.get_problem())
        geom=DropletGeometry(curv_radius=pr.R0,contact_angle=pr.contact_angle)
        self.mesh_mode="tris"
        self.default_resolution=0.1
        cl_factor=0.1
        p00=self.point(0,0)
        pr0=self.point(geom.base_radius,0,size=self.default_resolution*cl_factor)
        p0h=self.point(0,geom.apex_height)
       
        self.create_lines(p0h,"bubble_axis",p00,"bubble_substrate_Pt",pr0)
        self.circle_arc(pr0,p0h,through_point=(-geom.base_radius,0),name="interface")
        self.plane_surface("interface","bubble_substrate_Pt","bubble_axis",name="bubble")
       
        far_resolution=self.default_resolution*(pr.Lhost/pr.R0)
        pL0=self.point(pr.Lhost,0,size=far_resolution)
        p0L=self.point(0,pr.Lhost,size=far_resolution)
        pLL=self.point(pr.Lhost,pr.Lhost,size=far_resolution)
        pPG=self.point(pr.Pt_radius,0,size=self.default_resolution*cl_factor)
        pPGd = self.point(pr.Pt_radius, - pr.delta,size=self.default_resolution*cl_factor)
        pPGr = self.point(pr.Pt_radius + pr.delta,0,size=self.default_resolution*cl_factor)
        self.create_lines(pr0,"liquid_Pt",pPGr,"liquid_substrate",pL0,"liquid_side",pLL,"liquid_top",p0L,"liquid_axis",p0h)
        self.plane_surface("liquid_Pt","liquid_axis","liquid_substrate","liquid_side","interface","liquid_top",name="liquid")
 
        # add a substrate with thickness
        p0S=self.point(0,-pr.thickness)
        pLS=self.point(pr.Lhost,-pr.thickness,size=far_resolution)
        pPS=self.point(pr.Pt_radius,-pr.thickness)
        self.create_lines(p00,"substrate_axis",p0S,"substrate_base_Pt",pPS,"substrate_base",pLS,"substrate_side",pL0)
        self.create_lines(pPGr,"Pt_side_trangle",pPGd,"Pt_side",pPS)
        self.plane_surface("substrate_axis","substrate_base_Pt","Pt_side","liquid_Pt","bubble_substrate_Pt","Pt_side_trangle",name="Pt")
        self.plane_surface("Pt_side","substrate_base","substrate_side","liquid_substrate","Pt_side_trangle",name="substrate")
 
        # attach a remesher for meshing
        self.remesher = Remesher2d(self)
 
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
       
        self.create_lines(p0h,"bubble_axis",p00,"bubble_substrate_Pt",pr0)
        self.circle_arc(pr0,p0h,through_point=(-geom.base_radius,0),name="interface")
        self.plane_surface("interface","bubble_substrate_Pt","bubble_axis",name="bubble")
       
        far_resolution=self.default_resolution*(pr.Lhost/pr.R0)
        pL0=self.point(pr.Lhost,0,size=far_resolution)
        p0L=self.point(0,pr.Lhost,size=far_resolution)
        pLL=self.point(pr.Lhost,pr.Lhost,size=far_resolution)
        pPG=self.point(pr.Pt_radius,0,size=self.default_resolution/cl_factor)
        #pPGd = self.point(pr.Pt_radius, - pr.detal,size=self.default_resolution*cl_factor)
        #pPGr = self.point(pr.Pt_radius + pr.detal,0,size=self.default_resolution*cl_factor)
        self.create_lines(pr0,"liquid_Pt",pPG,"liquid_substrate",pL0,"liquid_side",pLL,"liquid_top",p0L,"liquid_axis",p0h)
        self.plane_surface("liquid_Pt","liquid_axis","liquid_substrate","liquid_side","interface","liquid_top",name="liquid")
 
        # add a substrate with thickness
        p0S=self.point(0,-pr.thickness)
        pLS=self.point(pr.Lhost,-pr.thickness,size=far_resolution)
        pPS=self.point(pr.Pt_radius,-pr.thickness)
        self.create_lines(p00,"substrate_axis",p0S,"substrate_base_Pt",pPS,"substrate_base",pLS,"substrate_side",pL0)
        self.create_lines(pPG,"Pt_side",pPS)
        self.plane_surface("substrate_axis","substrate_base_Pt","Pt_side","liquid_Pt","bubble_substrate_Pt",name="Pt")
        self.plane_surface("Pt_side","substrate_base","substrate_side","liquid_substrate",name="substrate")
 
        # attach a remesher for meshing
        self.remesher = Remesher2d(self)

class Current(Equations):
    def __init__(self, Jspace = "C2", Pspace = "C1"):
        super(Current, self).__init__()
        self.Jspace = Jspace
        self.Pspace = Pspace

    def define_fields(self):
        self.define_vector_field("J",self.Jspace) 
        self.define_scalar_field("phi",self.Pspace)

    def define_residuals(self):
        u,v = var_and_test("J")
        r,q = var_and_test("phi")
        self.add_residual(weak(u,v) + weak(r,div(v)*scale_factor("spatial")))
        self.add_residual(weak(div(u)*scale_factor("spatial"), q))

class CurrentNoFlux(InterfaceEquations):
     required_parent_type = Current # Must be attached to an interface of a Stokes equation

     def define_fields(self):
             self.define_scalar_field("noflux_lambda","C2")

     def define_residuals(self):
             # Binding variables
             l,ltest=var_and_test("noflux_lambda")
             u,utest=var_and_test("J")
             n=var("normal")
             self.add_residual(weak(dot(u,n),ltest)+weak(l,dot(utest,n)))

     # This will be called before the equations are numbered. This is the last chance to apply any pinning (i.e. Dirichlet conditions)
     def before_assigning_equations_postorder(self, mesh):
             self.pin_redundant_lagrange_multipliers(mesh, "noflux_lambda", "J")


# Problem class. Combining the equations, meshes, etc to a full problem
class SingleBubbleProblem(Problem):
    def __init__(self):
        super().__init__()
 
        # Define dimensions
        self.R0=0.65*milli*meter
        self.contact_angle=163*degree      
        self.Lhost=10*milli*meter
        self.thickness=2*milli*meter
        self.Pt_radius=1*milli*meter
        self.delta=0.1*milli*meter
 
        self.voltage_top = 0*volt
        self.voltage_base = -0.7*volt
        self.phi0 = -0.3*volt
        
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
 
        self.electric_conductivity = 40 *ampere/(volt*meter)
        # Contact line define
        #self.pinned_contact_line=True
        self.pinned_contact_line=False
       
        self.sliplength=0.01*micro*meter
       
        self.remesh_options=RemeshingOptions()
        self.remesh_options.active=True
 
        # Values at current time step
        #self.rho_air=var("gas_mass", domain="bubble")/var("gas_volume", domain="bubble") # density of air
        #self.mass_transfer_rate = self.diffusivity*dot(var("normal"),grad(var("c", domain="liquid"))) # mass transfer rate
 
    #def define_fields(self):
        #self.define_scalar_field("Q",space="C2")
    def decreaing_potential(self,r=var("coordinate_x")):
        phi_v = (1-r/self.Pt_radius)*self.voltage_base
        return phi_v
   
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
        self.set_scaling(temporal=1*milli*second)
        self.set_scaling(velocity=1*meter/second)
        self.set_scaling(pressure=self.sigma0/scale_factor("spatial"))
        self.set_scaling(c=self.c0)
        self.set_scaling(temperature=self.Troom)
        self.set_scaling(phiv=1*volt)
        self.set_scaling(thermal_equation=self.Troom*self.thermal_conductivity/self.R0**2)
        self.set_scaling(T=self.Troom)
        self.set_scaling(elec=1*ampere/(volt*meter))
        self.set_scaling(phiv=1*volt)
        self.set_scaling(Qscale=1*watt/meter**3)
        self.set_scaling(Jscale=1*ampere/meter**2)
        #self.set_scaling(J = 1 * ampere / meter**2)

        self.define_named_var(absolute_pressure=1*atm + var("pressure"))
 
        # add the mesh
        self+=SingleBubbleAxisymmMesh()
        self+=SingleBubblePlotter()
        
 
        geom=DropletGeometry(curv_radius=self.R0,contact_angle=self.contact_angle)
 
 
        #Mesh and output
        # Begin to assemble the equations for each of the three domains
        #bubble_eqs=MeshFileOutput()
        liquid_eqs=MeshFileOutput()
        Pt_eqs=MeshFileOutput()
        glass_eqs=MeshFileOutput()

        liquid_eqs+=LaplaceSmoothedMesh()
        Pt_eqs+=LaplaceSmoothedMesh()
        glass_eqs+=LaplaceSmoothedMesh()  

        #liquid_eqs+=Current()
        #liquid_eqs+=DirichletBC(phi=0)@"liquid_top"
        #liquid_eqs+=DirichletBC(phi=-1)@"liquid_Pt"
        #liquid_eqs+=CurrentNoFlux()@"interface"
        #liquid_eqs+=InitialCondition(J_x = 1, J_y=0)
        #liquid_eqs+=Scaling(J=scale_factor("Jscale"))+TestScaling(J=1/scale_factor("Jscale"))
        #liquid_eqs+=DirichletBC(J_y = 1, J_x = 0)@"liquid_Pt"
        #liquid_eqs+=DirichletBC(J_x = 0)@"liquid_axis"
        #liquid_eqs+=DirichletBC(J_x = 0)@"liquid_side"
        #liquid_eqs+=DirichletBC(J_y = 0)@"liquid_substrate"

        # Poisson (Laplace) equation for the electric potential:
        liquid_eqs+=PoissonEquation(name = "phiv",coefficient=scale_factor("spatial")**2)
        #liquid_eqs+=PoissonEquation(name = "phiv",coefficient=self.electric_conductivity*scale_factor("spatial")**2/scale_factor("elec"))
        liquid_eqs+=DirichletBC(phiv=self.voltage_top)@"liquid_top"
        #liquid_eqs+=NeumannBC(phiv = 0)@"liquid_side"
        #liquid_eqs+=NeumannBC(phiv = 0)@"liquid_substrate"
        #liquid_eqs+=NeumannBC(phiv = 0)@"interface"
        #liquid_eqs+=DirichletBC(phiv=self.voltage_base)@"liquid_Pt"
        #liquid_eqs+=DirichletBC(phiv=self.voltage_base)@"liquid_substrate"
        #liquid_eqs+=DirichletBC(phiv=self.decreaing_potential())@"liquid_Pt"
        liquid_eqs+=AxisymmetryBC()@"liquid_axis"
        #liquid_eqs+=InitialCondition(phi=-0.5)
        #liquid_eqs+=Scaling(phiv=scale_factor("phiv"))#+TestScaling(J=1/scale_factor("phiv"))
        liquid_eqs+=ConnectFieldsAtInterface("phiv")@"liquid_Pt"
        liquid_eqs+=ConnectMeshAtInterface()@"liquid_Pt"  
        #liquid_eqs+=DirichletBC(phiv=self.voltage_base)@"liquid_substrate/liquid_Pt"
        # remember to remove the _lagrange constrin when connecting to the bubble
        Pt_eqs+=PoissonEquation(name = "phiv",coefficient=1*self.electric_conductivity*scale_factor("spatial")**2/scale_factor("elec"))
        #Pt_eqs+=InitialCondition(phiv=self.phi0)
        Pt_eqs+=DirichletBC(phiv=self.voltage_base)@"substrate_base_Pt"
        Pt_eqs+=AxisymmetryBC()@"substrate_axis"
        #Pt_eqs+=NeumannBC(phiv = 0)@"bubble_substrate_Pt"
        #Pt_eqs+=NeumannBC(phiv = 0)@"Pt_side"

        #glass_eqs+=PoissonEquation(name = "phiv",coefficient=0.00000000001*self.electric_conductivity*scale_factor("spatial")**2/scale_factor("elec"))

        # Navier-Stokes equation in the liquid domain
        liquid_eqs+=NavierStokesEquations(mass_density=self.rho_liq,dynamic_viscosity=self.mu_liq,gravity=self.gravity)
        liquid_eqs+=AxisymmetryBC()@"liquid_axis"
        liquid_eqs+=LaplaceSmoothedMesh()
        liquid_eqs+=DirichletBC(velocity_x=0)@"liquid_side"
        liquid_eqs+=DirichletBC(mesh_y=0,velocity_y=0)@"liquid_substrate"
        liquid_eqs+=DirichletBC(mesh_y=0,velocity_y=0)@"liquid_Pt"
        liquid_eqs+=NavierStokesFreeSurface(surface_tension=self.surface_tension, static_interface=True)@"interface"  #, static_interface=True

        # In this case, we don't consider the expansion of the bubble

        # Source for temperature from the current

        
        #J=self.electric_conductivity*grad(var("phiv")) 
        #Q=dot(J,J)/self.electric_conductivity
        #liquid_eqs+=ProjectExpression(Q=dot(grad(var("phiv")),grad(var("phiv")))*self.electric_conductivity)
        #liquid_eqs+=Scaling(Q=scale_factor("Qscale"))+TestScaling(Q=1/scale_factor("Qscale"))
        #liquid_eqs+=ProjectExpression(J=self.electric_conductivity*grad(var("phiv")), field_type="vector" )
        #liquid_eqs+=Scaling(J=scale_factor("Jscale"))+TestScaling(J=1/scale_factor("Jscale"))
        #liquid_eqs+=DirichletBC(J_x = 0, J_y=0)@"liquid_Pt/liquid_substrate"
        #liquid_eqs+=ProjectExpression(Q=dot(var("J"),var("J"))/self.electric_conductivity)
        #liquid_eqs+=Scaling(Q=scale_factor("Qscale"))+TestScaling(Q=1/scale_factor("Qscale"))
        #Q=dot(var("J"),var("J"))/self.electric_conductivity
        #Q=dot(grad(var("phiv")),grad(var("phiv")))*self.electric_conductivity

        #liquid_eqs+=LocalExpressions(Q=dot(var("J"),var("J"))/self.electric_conductivity)
        liquid_eqs+=LocalExpressions(J=self.electric_conductivity*grad(var("phiv")) )
        #Pt_eqs+=LocalExpressions(J=self.electric_conductivity*grad(var("phiv")) )



        # Temperature and concentration advection-diffusion equation
        liquid_eqs+=AdvectionDiffusionEquations("temperature", diffusivity=self.thermal_diffusivity, space="C1", wind=var("velocity")) # , source = Q/ (self.mass_density*self.specific_heat_capacity)) #, source = Q/ (self.mass_density*self.specific_heat_capacity
        liquid_eqs+=AdvectionDiffusionInfinity(temperature=self.Troom)@"liquid_side"
        liquid_eqs+=AdvectionDiffusionInfinity(temperature=self.Troom)@"liquid_top"
        liquid_eqs+=InitialCondition(temperature=self.Troom)
        #liquid_eqs+=DirichletBC(temperature=self.Tsubstrate)@"liquid_substrate"
        liquid_eqs+=ConnectFieldsAtInterface("temperature")@"liquid_Pt"
        Pt_eqs+=AdvectionDiffusionEquations("temperature", diffusivity=self.thermal_diffusivity, space="C1", wind=0)
        Pt_eqs+=InitialCondition(temperature=self.Troom)



        #liquid_eqs+=LocalExpressions(j=grad(var("phiv")))
        #Pt_eqs+=LocalExpressions(j=grad(var("phiv")))

        self.add_equations(liquid_eqs@"liquid")
        self.add_equations(Pt_eqs@"Pt")
 
with SingleBubbleProblem() as p:
    #p.presolve_temperature()
    #p.run(10*second,outstep=True,temporal_error=1,startstep=0.001*second)#,start_step=0.001*second
    p.solve(max_newton_iterations=20)
    p.output()
    #p.output_at_increased_time()
    #p.run(10*second,outstep=True,temporal_error=1,start_step=0.001*second)
