#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
#Collision_V build data with v and no initial p. The variable is the same p, but it refers to v in the initial data
initial     = "Collision_V"

xrange = yrange = (-40., 40.)
zrange = (-40.0, 40.0)
Nx, Ny, Nz = 48, 48, 48

m0, q0, dm, dq = 0.5, 0.0, 10.0, 4.0
# dp  = (0.0, 0.0, 30.0)
dv  = (0.0, 0.0, 1.0)
x0  = (0.0, 0.0, 10.0)
sL  = (10., 10.)
sT  = (10., 10.)
mu  = dm^0.25
n   = 1

dt = 0.1 / mu
final_time = 20 / mu + dt

directory    = "/Users/mikel.sanchez/Documents/Chihuahua/data/collision_v/"
filename     = "symmetric-s-10-N-48-L-60-dv-1-dx-0.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "vx", "vz"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
