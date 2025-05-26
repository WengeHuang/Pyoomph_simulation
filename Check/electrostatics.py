from pyoomph import *
from pyoomph.expressions import *
from pyoomph.expressions.units import *
from pyoomph.expressions.phys_consts import epsilon_0 as epsilon_0
from pyoomph.equations.navier_stokes import StokesEquations

class ElectricPoissonEquation(Equations):
    def __init__(self,potential_name:str="V",charge_density:ExpressionOrNum=0,epsilon0:ExpressionOrNum=epsilon_0,epsilonr:ExpressionOrNum=1,epsilon:ExpressionOrNum=None,epsilon_reference:ExpressionOrNum=None):
        super().__init__()
        self.potential_name=potential_name
        self.charge_density=charge_density
        self.epsilon=epsilon
        if self.epsilon is None:
            self.epsilon=epsilon0*epsilonr
            if epsilon_reference is not None:
                self.epsilon_reference=epsilon_reference
            else:
                self.epsilon_reference=epsilon0
        else:
            raise RuntimeError("TODO: introduce some reference epsilon. If epsilon is a numerical constant, we could just use this. Otherwise, it likely will depend on whether it has a unit or not")
    
    def define_fields(self):
        self.define_scalar_field(self.potential_name,"C2",testscale=scale_factor("spatial")**2/(scale_factor(self.potential_name)*self.epsilon_reference))
        
    def define_residuals(self):
        V,Vtest=var_and_test(self.potential_name)
        self.add_weak(self.epsilon*grad(V),grad(Vtest)).add_weak(self.charge_density,Vtest)


class ElectricChangeTransportEquation(Equations):
    def __init__(self,conductivity:ExpressionOrNum,name:str="rho_e",wind:ExpressionNumOrNone=None,potential_name:Optional[str]=None):
        super().__init__()
        self.conductivity=conductivity
        self.name=name
        self.wind=wind
        self.potential_name=potential_name
        
    def define_fields(self):
        self.define_scalar_field(self.name,"C2",testscale=scale_factor("temporal")/scale_factor(self.name))
        
    def define_residuals(self):
        rhoe,rhoe_test=var_and_test(self.name)
        potname=self.potential_name
        if potname is None:
            potential_eqs=self.get_combined_equations().get_equation_of_type(ElectricPoissonEquation,always_as_list=True)
            if len(potential_eqs)!=1:
                raise RuntimeError("Cannot find a unique ElectricPoissonEquation to extract the name of the electric potential to use. Either set the potential_name or combine with a single ElectricPoissonEquation")
            potname=cast(ElectricPoissonEquation,potential_eqs[0]).potential_name
        u=self.wind
        if u is None:
            nseqs=self.get_combined_equations().get_equation_of_type(StokesEquations,always_as_list=True)
            if len(nseqs)==0:
                u=0
            elif len(nseqs)==1:
                u=var(cast(StokesEquations,nseqs[0]).velocity_name)
            else:
                raise RuntimeError("Cannot find a unique StokesEquations to extract the name of the velocity to use. Either set wind=... or combine with a single (Navier)StokesEquations")
        J=rhoe*u -self.conductivity*grad(var(potname))
        self.add_weak(partial_t(rhoe),rhoe_test).add_weak(-J,grad(rhoe_test))