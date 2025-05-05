from LN2_droplet import *

# Create the problem
problem=SingleBubbleLaserProblem()

# (optional) Attach a plotter to make nice images
problem+=SingleBubblePlotter(problem)

# Here, you can change some parameters

problem.R0=4*micro*meter # initial bubble radius (curvature radius)
problem.contact_angle=60*degree # contact angle of the bubble
problem.laser_power=12*milli*watt # total laser power
problem.laser_conversion_efficiency=0.4 # absorption coefficient
problem.laser_FWHM_width=10*micro*meter
problem.sliplength=0.01*micro*meter # slip length at the substrate
# Sticking coefficient of the mass transfer model 
# (0: no mass transfer, 1: maximum according to the kinetic gas theory, however, can also be higher)
# The lower it is, the more it will superheat/cool. The higher it is, the more the interface will stay at ~100°C (less Marangoni as consequence)
#problem.sticking_coefficient=0.1
problem.alpha.value=0.0

if False:
    lvl=4
    problem+=RefineToLevel(lvl)@"liquid/interface/liquid_substrate"
    problem+=RefineToLevel(lvl)@"bubble/interface/bubble_substrate"
    problem+=RefineToLevel(lvl)@"substrate/liquid_substrate/bubble_substrate"

T=1*micro*second # end time
outstep=0.1*micro*second # output interval

# presolve the temperature field
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