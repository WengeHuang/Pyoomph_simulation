from Leidenfrost import *

# Create the problem
problem=LeidenfrostProblem()

# (optional) Attach a plotter to make nice images
problem+=LeidenfrostPlotter(problem)

# Here, you can change some parameters

problem.droplet_height=1.1*milli*meter
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

T=10*milli*second # end time
outstep=0.1*milli*second # output interval

# presolve the temperature field
#problem.solve(max_newton_iterations=20)
problem.presolve_temperature()
# and run the simulation with temporal adaptivity
problem.run(T,outstep=outstep,temporal_error=1)


#problem.presolve_temperature()
#problem.run(100*micro*second, outstep = False, startstep = 0.1*micro*second, temporal_error=1)
#problem.solve()
#problem.output()

#while True:
#    problem.go_to_param(alpha=problem.alpha+1)
#    problem.remesh_handler_during_continuation()
#    problem.output_at_increased_time()

#ds = 1
#while problem.alpha.value<50:
#    ds = problem.arclength_continuation("alpha", ds)
#    problem.remesh_handler_during_continuation()
#    problem.output_at_increased_time()