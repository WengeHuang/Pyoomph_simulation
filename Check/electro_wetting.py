from electrostatics import *
from pyoomph.equations.navier_stokes import *
from pyoomph.utils.dropgeom import DropletGeometry
from pyoomph.equations.ALE import *


class DropletMesh(GmshTemplate):
    def define_geometry(self):
        # Just a simple mesh of a droplet on a substrate with gas around it
        self.default_resolution=0.05
        pr=cast(ElectroWettingProblem,self.get_problem())
        p0=self.point(0,0)
        geom=DropletGeometry(volume=pr.volume,contact_angle=pr.theta0)
        pr0=self.point(geom.base_radius,0,size=-0.3)
        p0h=self.point(0,geom.apex_height)
        self.circle_arc(pr0,p0h,through_point=(-geom.base_radius,0),name="droplet_gas")
        self.line(p0,pr0,name="droplet_substrate")
        self.line(p0,p0h,name="droplet_axis")
        self.plane_surface("droplet_gas","droplet_axis","droplet_substrate", name="droplet")
        pR0=self.point(pr.Rfar_factor*geom.base_radius,0,size=-pr.Rfar_factor)
        self.plane_surface("droplet_substrate",*self.create_lines(pr0,"gas_substrate",pR0,"substrate_far",self.point(pr.Rfar_factor*geom.base_radius,-pr.Hsubstrate,size=-pr.Rfar_factor),"substrate_bottom",self.point(0,-pr.Hsubstrate),"substrate_axis",p0),name="substrate")
        p0R=self.point(0,pr.Rfar_factor*geom.base_radius,size=-pr.Rfar_factor)
        self.plane_surface("droplet_gas","gas_substrate",self.circle_arc(pR0,p0R,center=p0,name="gas_far"),self.line(p0R,p0h,name="gas_axis"),name="gas")
                                

class ElectroWettingProblem(Problem):
    def __init__(self):
        super().__init__()
        self.volume=1*micro*liter
        self.theta0=90*degree        
        self.V=self.define_global_parameter(V=0.1)
        self.voltage=self.V*volt
        self.Rfar_factor=4
        self.Hsubstrate=0.1*milli*meter
        self.rho=1000*kilogram/meter**3
        self.mu=1*milli*pascal*second
        self.sigma=72.8*milli*newton/meter
        self.lslip=100*micro*meter
        # Relative permittivities        
        self.eps_r_gas=1 # air 
        self.eps_r_drop=80  # water DC 20°C
        self.eps_r_substrate=3 # Acrylic
        
        
    def define_problem(self):
        geom=DropletGeometry(volume=self.volume,contact_angle=self.theta0)
        
        self.set_scaling(spatial=square_root(self.volume,3),V=volt,temporal=1*second,velocity=scale_factor("spatial")/scale_factor("temporal"),pressure=2*self.sigma/geom.curv_radius)
        self.set_coordinate_system("axisymmetric")
        self+=DropletMesh()
        
        # Droplet bulk equations, Navier-Stokes on moving mesh with Poisson equation for the electric potential
        deqs=MeshFileOutput()+NavierStokesEquations(mass_density=self.rho,dynamic_viscosity=self.mu)
        deqs+=ElectricPoissonEquation(epsilonr=self.eps_r_drop,epsilon_reference=epsilon_0)
        deqs+=PseudoElasticMesh()+EnforceVolumeByPressure(self.volume)
        
        # BCs for the droplet
        deqs+=AxisymmetryBC()@"droplet_axis"
        deqs+=(DirichletBC(mesh_y=0,velocity_y=0)+NavierStokesSlipLength(self.lslip))@"droplet_substrate"        
        gradV=grad(var("V",domain="droplet")) # Bulk gradient at the interface                        
        electric_surface_force=self.eps_r_drop*epsilon_0/2*dot(gradV,gradV)
        deqs+=TextFileOutput()@"droplet_gas"
        # Add the electric surface force 
        deqs+=WeakContribution(electric_surface_force,-dot(var("normal"),testfunction("velocity")))@"droplet_gas"
        
        # Bulk effects of the electric field?? No idea...
        #deqs+=WeakContribution(self.sigma,div(testfunction("velocity")))@"droplet_substrate"
        #E=-grad(var("V"))
        #TE=self.eps_r_drop*epsilon_0*(dyadic(E,E)-1/2*dot(E,E)*identity_matrix())
        #deqs+=WeakContribution(TE,grad(testfunction("velocity")))
        
        # Fix the voltage at the free interface. We use and EnforcedDirichlet condition, which is better than a DirichletBC if you use continuation
        # Also impose the contact angle condition without voltage. This balances the forces at the contact line like a Neumann triangle.
        deqs+=(EnforcedDirichlet(V=self.voltage)+NavierStokesFreeSurface(surface_tension=self.sigma)+NavierStokesContactAngle(contact_angle=self.theta0)@"droplet_substrate")@"droplet_gas"                
        
        # Substrate equations, Poisson equation for the electric potential and moving mesh
        seqs=MeshFileOutput()+PseudoElasticMesh()    
        seqs+=ElectricPoissonEquation(epsilonr=self.eps_r_substrate,epsilon_reference=epsilon_0)
        # BCs for the substrate
        seqs+=AxisymmetryBC()@"substrate_axis"
        seqs+=DirichletBC(mesh_y=0)@"droplet_substrate"
        seqs+=DirichletBC(mesh_y=0)@"gas_substrate"
        seqs+=DirichletBC(mesh_y=-self.Hsubstrate)@"substrate_bottom"
        seqs+=DirichletBC(mesh_x=self.Rfar_factor*geom.base_radius)@"substrate_far"        
        seqs+=DirichletBC(V=0)@"substrate_bottom"
        
        # Gas equations, Poisson equation for the electric potential and moving mesh
        geqs=MeshFileOutput()+PseudoElasticMesh()    
        geqs+=ElectricPoissonEquation(epsilonr=self.eps_r_gas,epsilon_reference=epsilon_0)
        geqs+=AxisymmetryBC()@"gas_axis"
        geqs+=PinMeshCoordinates()@"gas_far"
        geqs+=DirichletBC(mesh_y=0)@"gas_substrate"
        
        # Connections
        deqs+=(ConnectMeshAtInterface()+ConnectFieldsAtInterface("V"))@"droplet_substrate"
        deqs+=(ConnectMeshAtInterface()+ConnectFieldsAtInterface("V"))@"droplet_gas"
        geqs+=(ConnectMeshAtInterface()+ConnectFieldsAtInterface("V"))@"gas_substrate"
        self+=deqs@"droplet"+seqs@"substrate"+geqs@"gas"
        
        # Make some ODE to measure the contact angle
        theta_deg,theta_deg_test=self.add_global_dof("theta_deg",initial_condition=self.theta0/degree)
        self+=WeakContribution(theta_deg*degree-acos(var("normal_x")),theta_deg_test)@"droplet/droplet_gas/droplet_substrate"
        theta_deg_lipp=acos(cos(self.theta0)+self.eps_r_substrate*epsilon_0*self.voltage**2/(2*self.sigma*self.Hsubstrate))
        self.add_global_dof("theta_deg_lippmann",initial_condition=self.theta0/degree,equation_contribution=var("theta_deg_lippmann")*degree-theta_deg_lipp)
        self+=ODEFileOutput("contact_angle.txt")@"globals"
        
        

pr=ElectroWettingProblem()
if False:
    pr.max_refinement_level=4
    pr+=RefineToLevel()@"droplet/droplet_gas/droplet_substrate"
    pr+=RefineToLevel()@"gas/droplet_gas/gas_substrate"
    pr+=RefineToLevel()@"substrate/droplet_substrate/droplet_gas"
    pr+=RefineToLevel()@"substrate/gas_substrate/droplet_gas"
pr.solve(spatial_adapt=pr.max_refinement_level)
pr.output()
pr.go_to_param(V=100)
pr.output_at_increased_time()
pr.go_to_param(V=200)
pr.output_at_increased_time()