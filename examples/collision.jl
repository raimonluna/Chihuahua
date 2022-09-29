#########################
# PARAMETERS
#########################

xrange = yrange = zrange = (-25.0, 25.0)
Nx = Ny = Nz = 64

m0, q0, dm, dq = 0.5, 0.0, 10.0, 6.0 
dp = (0.0, 0.0, 30.0)
x0 = (2.0, 0.0, 10.0)
s  = 10.0
mu = dm^0.25

dt = 0.1 / mu
final_time = 20 / mu + dt

out_file     = pwd() * "/collision_output_2.m"
out_funcs_xz = ["mass", "charge"]
out_every    = 10

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
