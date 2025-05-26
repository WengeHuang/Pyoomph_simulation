from Leidenfrost import *

# Create the problem
problem=LeidenfrostProblem()

# (optional) Attach a plotter to make nice images
problem+=LeidenfrostPlotter(problem)

# Here, you can change some parameters

problem.droplet_height=1.01*milli*meter
problem.droplet_radius=1*milli*meter
# Sticking coefficient of the mass transfer model 
# (0: no mass transfer, 1: maximum according to the kinetic gas theory, however, can also be higher)
# The lower it is, the more it will superheat/cool. The higher it is, the more the interface will stay at ~100°C (less Marangoni as consequence)
#problem.sticking_coefficient=0.1


if False:
    lvl=4
    problem+=RefineToLevel(lvl)@"liquid/interface/liquid_substrate"
    problem+=RefineToLevel(lvl)@"bubble/interface/bubble_substrate"
    problem+=RefineToLevel(lvl)@"substrate/liquid_substrate/bubble_substrate"

T=3000*milli*second # end time
outstep=1*milli*second # output interval

# presolve the temperature field
#problem.solve(max_newton_iterations=20)
problem.presolve_temperature()
#problem.solve()
#problem.output()
# and run the simulation with temporal adaptivity
problem.run(T,outstep=outstep,temporal_error=1,out_initially=False)# , out_initially=False
