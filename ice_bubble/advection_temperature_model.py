import sys
sys.path.append('../..')
from problem_def.problem import *

# Create a plotter object
class Plotter(MatplotlibPlotter):
    def __init__(self, problem: Problem, filetrunk: str = "plot_{:05d}", fileext: str | List[str] = "pdf", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

    def define_plot(self):
        pr=cast("BubbleFreezingProblem",self.get_problem())
        # h1=pr.get_ode("globals").get_value("ice_front_pos_y")
        h1=0
        h2=pr.channel_height*3/4
        r=pr.channel_width/1.5
        self.dpi=300

        # Set view
        self.set_view(xmin=0,xmax=r,ymin=h1,ymax=1*(h1+h2))

        # Colorbars
        cbT=self.add_colorbar("c [kg/m$^3$]", position="top right", cmap="coolwarm")
        cbT.invisible=True
        cb_streams = self.add_colorbar("stream_sign", cmap="Greys", position="center")  # black to white
        cb_streams.Ndisc = 2  # but only one discrete colors (black)
        cb_streams.symmetrize_min_max = True  # Symmetrize the colorbar
        cb_streams.invisible = True  # Do not show the colorbar
        streams = self.add_plot("liquid/velocity", mode="streamlines", linewidths=0.5, linecolor="white")

        # Add plots
        self.add_plot("liquid/c", colorbar=cbT, transform=[None,"mirror_x"])
        self.add_plot("liquid/ice_front",transform=[None,"mirror_x"],linecolor="purple")
        self.add_plot("liquid/bubble",transform=[None,"mirror_x"],linecolor="green")
        self.add_plot("liquid/T", colorbar=cbT, transform="mirror_x")   ##comment out 
        
        # Add text
        # t=round(float(pr.get_current_time()/second),3)
        # radius=round(float(pr.get_ode("globals").get_value("bubble_radius")/pr.bubble_radius),3)
        # txt = self.add_text("$R$ = {}".format(radius), "top center", bbox=dict(facecolor='wheat', alpha=1, boxstyle="round"))
        # txt.ymargin-=0.0254

with BubbleFreezingProblem() as p:
    p.plotter = Plotter(p)
    p.pure_diffusion_static_uniform_c_model = False
    p.pure_diffusion_static_model = False
    p.pure_diffusion_moving_model = False
    p.advection_temperature_model = False     ## comment out 
    p.advection_fake_temperature_model = True   
    p.presolve_gas_phase()
    p.run(endtime=10*second, startstep=0.01*second, outstep=0.01*second, temporal_error=1)
    p.output_at_increased_time()