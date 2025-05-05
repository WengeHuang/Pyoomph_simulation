from pyoomph import *
from pyoomph.generic.problem import Problem
from pyoomph.meshes.meshdatacache import MeshDataEigenModes
from pyoomph.output.plotting import *
from pyoomph.typings import List, Optional, Union
from pyoomph.expressions.units import *

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


# Create a plotter object
class ZoomedPlotter(MatplotlibPlotter):
    def __init__(self, problem: Problem, filetrunk: str = "zoomed_{:05d}", fileext: str | List[str] = "png", eigenvector: int | None = None, eigenmode: MeshDataEigenModes = "abs", add_eigen_to_mesh_positions: bool = True, position_eigen_scale: float = 1):
        super().__init__(problem, filetrunk, fileext, eigenvector, eigenmode, add_eigen_to_mesh_positions, position_eigen_scale)

    def define_plot(self):
        pr=cast("BubbleFreezingProblem",self.get_problem())
        y_offset=pr.get_ode("globals").get_value("bubble_y")
        bubble_radius=pr.get_ode("globals").get_value("bubble_radius")

        # Set view
        xmin=-1.5*bubble_radius
        xmax=1.5*bubble_radius
        ymin=-1.5*bubble_radius+y_offset
        ymax=1.5*bubble_radius+y_offset
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
        if pr.advection_temperature_model or pr.advection_fake_temperature_model:
            self.add_plot("liquid/T", colorbar=cbT, transform="mirror_x")

        # Add time label
        self.add_time_label("top center")