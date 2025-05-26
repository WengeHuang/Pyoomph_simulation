from pyoomph import *
from pyoomph.expressions import *
from pyoomph.expressions.units import *
# Load tools for periodic driving response and text file output
from pyoomph.utils.periodic_driving_response import *
from pyoomph.utils.num_text_out import *

# Driven damped harmonic oscillator
class DampedHarmonicOscillatorEquations(ODEEquations):
    def __init__(self,omega0,delta,driving):
        super().__init__()
        self.omega0,self.delta,self.driving=omega0,delta,driving

    def define_fields(self):
        # Must be formulated first order here
        #self.define_ode_variable("y",testscale=scale_factor("temporal")**2/scale_factor("y"))
        #self.define_ode_variable("ydot",scale=scale_factor("y")/scale_factor("temporal"),testscale=scale_factor("temporal")/scale_factor("y"))

        self.define_ode_variable("y",testscale=scale_factor("temporal")/scale_factor("y") )
        self.define_ode_variable("ydot",scale=scale_factor("y")/scale_factor("temporal"),testscale=scale_factor("temporal")**2/scale_factor("y"))


    def define_residuals(self):
        y,y_test=var_and_test("y")
        ydot,ydot_test=var_and_test("ydot")
        #self.add_weak(partial_t(y)-ydot,ydot_test)
        #self.add_weak(partial_t(ydot)+self.delta*ydot +self.omega0**2*y-self.driving,y_test)
        self.add_weak(partial_t(y)-ydot,y_test)
        self.add_weak(partial_t(ydot)+self.delta*ydot +self.omega0**2*y-self.driving,ydot_test)


class DampedHarmonicOscillatorProblem(Problem):
    def __init__(self):
        super().__init__()
        self.omega0=1/second
        self.delta=0.1/second
        # Default driving
        self.driving=meter/second**2 *cos(0.3/second*var("time"))

    def define_problem(self):
        self.set_scaling(y=2*meter,temporal=1*second)
        eqs=DampedHarmonicOscillatorEquations(self.omega0,self.delta,self.driving)
        eqs+=ODEFileOutput()
        self+=eqs@"oscillator"

with DampedHarmonicOscillatorProblem() as problem:
     # Trivial, but long way: integrate in time, extract response manually from the output
     problem.run(500*second,outstep=0.1*second)