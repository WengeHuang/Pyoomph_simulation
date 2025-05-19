from pyoomph import *
from pyoomph.expressions import *  # Import grad & weak
from pyoomph.equations.poisson import *
from pyoomph.output.plotting import MatplotlibPlotter # Plotting tools

 
# Plotter calss. Will generate figures
class MixPoissonPlotter(MatplotlibPlotter):
    def define_plot(self):
        self.background_color = "darkgrey"
        self.set_view(-1, 0, 0, 1)  # x_min, ymin, x_max, y_max of the view window
        cb_c=self.add_colorbar("phi/J",cmap="Blues",position="top left")
        self.add_plot("domain/J",colorbar=cb_c,transform="mirror_x")
        self.add_plot("domain/J",mode="streamlines",transform="mirror_x")







class MixedPoissonEquation(Equations):
    def __init__(self, *, name="u",  source=0):
        super(MixedPoissonEquation, self).__init__()
        self.source=source
        

    def define_fields(self):
        self.define_scalar_field("phi", "D1")
        self.define_vector_field("J","C2")

    def define_residuals(self):
        phi,phitest=var_and_test("phi")
        J,Jtest=var_and_test("J")
        a=weak(J,Jtest)+weak(phi,div(Jtest))+weak(div(J),phitest)+weak(self.source,phitest)
        self.add_residual(a)
        
class MixedPoissonTest(Problem):
    def __init__(self):
        super().__init__()
        self.source=0
        self.phi0=1
        x,y=var(["coordinate_x","coordinate_y"])
        self.g=sin(5*x)
        self.f=10*exp(-((x-0.5)**2+(y-0.5)**2)/0.02)
        self.use_mixed=False
        
        
    def define_problem(self):
        self+=RectangularQuadMesh(N=40)
        self+=MixPoissonPlotter()
        eqs=MeshFileOutput()
        eqs+=MeshFileOutput("dg_output",discontinuous=True)            
        if not self.use_mixed:            
            eqs+=PoissonEquation(name="phi",source=self.f)
            eqs+=ProjectExpression(J=grad(var("phi")),field_type="vector",space="D1")
            eqs+=NeumannBC(phi=-self.g)@"top"
            eqs+=NeumannBC(phi=-self.g)@"bottom"
            eqs+=DirichletBC(phi=self.phi0)@"left"
            eqs+=DirichletBC(phi=self.phi0)@"right"
        else:
            eqs+=MixedPoissonEquation(source=self.f)                                   
            # Opposed to the normal Poisson equation implementation, we now have to set DirichletBCs for fluxes J (which are usually Neumann terms)
            # And NeumannBCs to prescibe phi (which are usually Dirichlet terms)
            eqs+=DirichletBC(J_y=self.g)@"top"
            eqs+=DirichletBC(J_y=-self.g)@"bottom"
            eqs+=NeumannBC(J=-var("normal")*self.phi0)@"left"
            eqs+=NeumannBC(J=-var("normal")*self.phi0)@"right"
            
            
        self+=eqs@"domain"
        
with MixedPoissonTest() as problem:
    #problem.set_coordinate_system(axisymmetric)
    problem.solve()
    problem.output()
        
