#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
#Collision_V build data with v and no initial p. The variable is the same p, but it refers to v in the initial data
initial     = "Collision"

xrange = yrange = (-40., 40.)
zrange = (-40.0, 40.0)
Nx, Ny, Nz = 40, 40, 120

m0, q0, dm, dq = 0.5, 0.0, 10.0, 3.0
dp  = (0.0, 0.0, 30.0)
x0  = (6.0, 0.0, 10.0)
sL  = (10., 10.)
sT  = (10., 10.)
mu  = dm^0.25
n   = 1

# dt = 0.1 / mu
dt  = 0.05
# final_time = 15 / mu + dt
final_time = 40.

directory    = "/Users/mikel.sanchez/Dropbox/Projects/data_Chihuahua/collisions/symmetric_shocks/q-3/"
filename     = "Nxy-40-Nz-120-dx-6-dp-30-dt-0.05.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "vx", "vy", "vz"]
# out_every    = Int((dt * mu)^-1)
out_every    = Int(dt^-1/2)
# out_every    = 1


#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
