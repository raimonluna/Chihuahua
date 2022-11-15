#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
#Collision_V build data with v and no initial p. The variable is the same p, but it refers to v in the initial data
initial     = "Single"

xrange = yrange = (-50., 50.)
zrange = (-50., 50.)
Nx, Ny, Nz = 50, 50, 150

m0, q0, dm, dq = 0.5, 0.0, 10.0, 3.0
dp  = (0.0, 0.0, 30.0)
x0  = (0.0, 0.0, 10.0)
sL  = 10.
sT  = 50.
mu  = dm^0.25

# dt = 0.1 / mu
dt  = 0.05
# final_time = 15 / mu + dt
final_time = 30.

directory    = "/Users/mikel.sanchez/Dropbox/Projects/data_Chihuahua/single_gaussian/"
filename     = "Nxy-50-Nz-150-dp-30-dt-0.05.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "vx", "vy", "vz"]
# out_every    = Int((dt * mu)^-1)
out_every    = Int(dt^-1/2)
# out_every    = 1


#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
