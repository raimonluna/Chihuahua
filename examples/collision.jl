#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
initial     = "Collision"

xrange = yrange = (-30., 30.)
zrange = (-30.0, 30.0)
Nx, Ny, Nz = 30, 32, 48

m0, q0, dm, dq = 0.5, 0.0, 50.0, 15.0
dp  = (0.0, 0.0, 30.0)
x0  = (0.0, 0.0, 4.0)
sL  = (4., 4.)
sT  = (4., 4.)
mu  = dm^0.25
n   = 1

dt = 0.1 / mu
final_time = 20 / mu + dt

directory    = "/Users/mikel.sanchez/Documents/Chihuahua/data/"
filename     = "a.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "pz", "vz"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
