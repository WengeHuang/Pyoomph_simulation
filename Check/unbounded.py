# Load pyoomph and some equations
from pyoomph import *
from pyoomph.equations.navier_stokes import *
from pyoomph.equations.poisson import *
from pyoomph.equations.ALE import *

from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools

 
# Plotter calss. Will generate figures
class DropletDeformPlotter(MatplotlibPlotter):
    def define_plot(self):
        # Get the problem object
        pr=cast(DropDeformProblem,self.get_problem())
        xrange = pr.Ro # view range
        yrange = pr.Ro
        self.background_color = "darkgrey"
        self.set_view(-2*xrange, -2*yrange, 2*xrange, 2*xrange)  # -x_min, -ymin, x_max, y_max of the view window
        #self.set_view(-0.5*xrange, -0.05*xrange, 0, 0.4*xrange)
        #cb_T = self.add_colorbar("temperature [°C]", offset=-273.15, position = "bottom right")
        cb_u = self.add_colorbar("J", position = "bottom left", cmap = "viridis")

        #self.add_plot("liquid/temperature",colorbar=cb_T,transform=[None,"mirror_x"])
       
        self.add_plot("droplet/velocity",colorbar=cb_u,transform=[None,"mirror_x"])
        self.add_plot("outer/velocity",colorbar=cb_u,transform=[None,"mirror_x"])

        self.add_plot("droplet/velocity",mode="streamlines",transform=[None,"mirror_x"])
        self.add_plot("outer/velocity",mode="streamlines",transform=[None,"mirror_x"])
        self.add_plot("outer/interface",transform=[None,"mirror_x"])
        CaE=self.add_text("CaE={:5g}".format(pr.CaE.value),position="bottom right",bbox=dict(boxstyle='round', facecolor='wheat', alpha=1),textsize=16)
        tl=self.add_time_label("bottom center")

# Just half of the domain due to symmetry
class DropDeformMesh(GmshTemplate):    
    def define_geometry(self):                
        pr=cast(DropDeformProblem,self.get_problem())
        W,H=pr.W,pr.H
        sizeD=0.05
        self.default_resolution=maximum(0.05*H,sizeD)
        self.mesh_mode="tris"
        # Points
        pN=self.point(0,pr.Ro,size=sizeD) # North pole of the droplet
        pS=self.point(0,-pr.Ro,size=sizeD) # South pole of the droplet
        pE=self.point(pr.Ro,0,size=sizeD) # East side pole of the droplet
        # Lines of the droplet
        self.circle_arc(pN,pE,through_point=pS,name="interface")
        self.circle_arc(pS,pE,through_point=pN,name="interface")
        #self.circle_arc(pN, pS, through_point=pE, name = "interface")
        # Here we need to define the interface twice since 
        # since it only creates the arc from the starting point to the ending point
        self.line(pN,pS,name="droplet_symm")
        # Create a line loop for the outer walls of the box
        walls=self.create_lines(pS,"outer_symm",(0,-H/2),"bottom",(W/2,-H/2),"right",(W/2,H/2),"top",(0,H/2),"outer_symm",pN)        
        # Surfaces
        self.plane_surface(*walls,"interface",name="outer")
        #The * operator unpacks a list or tuple into individual arguments when calling a function.
        self.plane_surface("interface","droplet_symm",name="droplet")



# Just a storage class for the domain information
class DomainInformation:
    def __init__(self,k,eps,mu) -> None:
        self.k,self.eps,self.mu=k,eps,mu

class DropDeformProblem(Problem):
    def __init__(self):
        super().__init__()
        # Geometric parameters
        self.H=20
        self.W=20
        self.Ro=1
        
        # Physical parameters 
        self.k1=1
        self.rho1=0.001
        self.mu1=1
        self.eps1=1
        self.gamma=0.1        
        self.CaE=0.5
        
        # Ratios
        self.R=1.75
        self.Q=3.5
        self.beta=1
        
        
    def get_deformation(self):
        # Calculate the deformation D
        imesh=self.get_cached_mesh_data("droplet/interface")
        segs,_=imesh.get_interface_line_segments()
        seg=segs[0]
        x,y=imesh.get_coordinates()
        b=abs(y[seg[0]]-y[seg[-1]])        
        a=2*numpy.amax(x) # TODO: Better way to get the radius
        import scipy.interpolate
        if y[seg[0]]>y[seg[-1]]:
            seg=seg[::-1]        
        inter=scipy.interpolate.UnivariateSpline(y[seg],x[seg],k=4,s=0)        
        y0=inter.derivative().roots()[0]
        a=2*inter(y0)
        return (b-a)/(b+a)
        
        
    def get_deformations(self):
        # Calculate the deformation D according to different approximations
        Dnum=self.get_deformation()
        Dtaylor=9/16*self.CaE/(2+self.R)**2*(1+self.R**2-2*self.Q+3/5*(self.R-self.Q)*(2+3*self.beta)/(1+self.beta))
        Dfeng=(self.R**2+self.R+1-3*self.Q)/(3*(1+self.R)**2)*self.CaE
        fd=(self.R**2+3/2*self.R+1-7/2*self.Q)/self.R**2
        K1=9/16*(self.R**2/(2+self.R)**2)*fd
        K2=((139*self.R-154)*self.R**2*fd + 92*(self.R**2+2*self.R-(self.R+2)*self.Q) )/(80*(2+self.R)**3)*K1
        D2nd=K1*self.CaE+K2*self.CaE**2
        return float(Dnum),float(Dtaylor),float(Dfeng),float(D2nd)
        
        
    def define_problem(self):
        # Define the problem
        self+=DropDeformMesh()
        self+=DropletDeformPlotter()
        droplet=DomainInformation(self.k1,self.eps1,self.mu1)
        outer=DomainInformation(self.k1/self.R,self.eps1/self.Q,self.mu1/self.beta)
        
        for domain,info in zip(["droplet","outer"],[droplet,outer]):
            eqs=MeshFileOutput()
            eqs+=PoissonEquation("phi",coefficient=info.k)
            #eqs+=NavierStokesEquations(mass_density=self.rho1,dynamic_viscosity=info.mu)
            eqs+=StokesEquations(dynamic_viscosity=info.mu)
            eqs+=LaplaceSmoothedMesh()
            self+=eqs@domain
            
        Einfty=square_root(self.CaE*self.gamma/(self.Ro*self.eps1/self.Q))
        self+=DirichletBC(phi=-self.H/2*Einfty,mesh_x=True,mesh_y=-self.H/2,velocity_x=0,velocity_y=0)@"outer/bottom"
        self+=DirichletBC(phi=self.H/2*Einfty,mesh_x=True,mesh_y=self.H/2,velocity_x=0,velocity_y=0)@"outer/top"
        self+=DirichletBC(mesh_y=True,mesh_x=self.W/2,velocity_x=0,velocity_y=0)@"outer/right"
        self+=DirichletBC(mesh_x=0,velocity_x=0)@"outer/outer_symm"        
        self+=DirichletBC(mesh_x=0,velocity_x=0)@"droplet/droplet_symm"        
        
        Edrop=-grad(var("phi",domain="droplet"))
        Eouter=-grad(var("phi",domain="outer"))
        Tdrop=self.eps1*dyadic(Edrop,Edrop)-1/2*self.eps1*dot(Edrop,Edrop)*identity_matrix()
        Touter=self.eps1/self.Q*dyadic(Eouter,Eouter)-1/2*self.eps1/self.Q*dot(Eouter,Eouter)*identity_matrix()
        
        
        ieqs=ConnectFieldsAtInterface(["phi","velocity_x","velocity_y"])+ConnectMeshAtInterface()
        ieqs+=NavierStokesFreeSurface(surface_tension=self.gamma)
        ieqs+=WeakContribution(matproduct(Tdrop-Touter,var("normal")),testfunction("velocity"))
        self+=ieqs@"droplet/interface"
        
        avgp,avgp_test=self.add_global_dof("avgp")
        self+=(WeakContribution(var("pressure"),avgp_test)+WeakContribution(avgp,testfunction("pressure")))@"outer"
        self+=(WeakContribution(var("pressure"),avgp_test)+WeakContribution(avgp,testfunction("pressure")))@"droplet"
        
        avgy,avgy_test=self.add_global_dof("avgy")
        #self+=(WeakContribution(var("coordinate_y"),avgy_test)+WeakContribution(avgy,testfunction("velocity_y")))@"droplet"
        self+=(WeakContribution(var("mesh_y"),avgy_test)+WeakContribution(avgy,testfunction("velocity_y")))@"droplet"

        fixV,fixV_test=self.add_global_dof("fixV",equation_contribution=-pi/2*self.Ro)
        self+=(WeakContribution(fixV,testfunction("pressure"))+WeakContribution(1,fixV_test))@"droplet"
        #self+=(WeakContribution(fixV,testfunction("pressure")))@"droplet"
        
        
            
            
with DropDeformProblem() as problem:    
    problem.CaE=problem.define_global_parameter(CaE=0.5)
    
    problem.solve()    
    problem.go_to_param(CaE=0.01)
    from pyoomph.utils.num_text_out import NumericalTextOutputFile
    outf=NumericalTextOutputFile(problem.get_output_directory("deformations.txt"),header=["CaE","Dnum","Dtaylor","Dfeng","D2nd"])
    for CaE in numpy.linspace(0.01,1,100):
        problem.go_to_param(CaE=CaE)
        problem.output_at_increased_time()
        outf.add_row(CaE,*problem.get_deformations())
    problem.output()

        