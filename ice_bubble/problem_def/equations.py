from pyoomph import *
from pyoomph.expressions import *
from pyoomph.expressions.units import *

class IceFrontSpeed(InterfaceEquations):
    def __init__(self, ice_front_speed: float = 1e-5*meter/second):
        super(IceFrontSpeed, self).__init__()
        self.ice_front_speed=ice_front_speed

    def define_fields(self):
        self.define_scalar_field("_lagr_interf_speed","C2",scale=1/test_scale_factor("mesh"),testscale=scale_factor("temporal")/scale_factor("spatial"))

    def define_residuals(self):
        speed=self.ice_front_speed
        n=var("normal")
        x,xtest=var_and_test("mesh")
        l,ltest=var_and_test("_lagr_interf_speed")
        self.add_residual(weak(dot(partial_t(x),n)-speed,ltest))
        self.add_residual(weak(l,dot(xtest,n)))