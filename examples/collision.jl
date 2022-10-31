#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
#Collision_V build data with v and no initial p. The variable is the same p, but it refers to v in the initial data
initial     = "Collision"

xrange = yrange = (-30., 30.)
zrange = (-30.0, 30.0)
Nx, Ny, Nz = 48, 48, 48

m0, q0, dm, dq = 0.2, 0.0, 50.0, 15.0
dp  = (0.0, 0.0, 40.0)
x0  = (0.0, 0.0, 5.0)
sL  = (4., 4.)
sT  = (4., 4.)
mu  = dm^0.25
n   = 1

dt = 0.1 / mu
final_time = 20 / mu + dt

directory    = "/Users/mikel.sanchez/Documents/Chihuahua/data/collision/"
filename     = "symmetric-s-4-N-48-L-40-dp-40-dx-0.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "vx", "vz"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
