#########################
# PARAMETERS
#For opposite charge choose n to -1
#########################
derivatives = "spectral"
#Collision_V build data with v and no initial p. The variable is the same p, but it refers to v in the initial data
initial     = "Collision_1D"

xrange = yrange = (-0.5, 0.5)
zrange = (-40.0, 40.0)
Nx, Ny, Nz = 4, 4, 150

m0, q0, dm, dq = 0.5, 0.0, 10.0, 6.0
dp             = 30.0
dz             = 20.0
s              = 10.
mu             = dm^0.25

# dt = 0.1 / mu
dt = 0.025
# final_time = 15 / mu + dt
final_time = 40

directory    = "/Users/mikel.sanchez/Dropbox/Projects/data_Chihuahua/noise_tests/1D/"
filename     = "1D_N-150-dp-30-dq-6-dt-0.025-match-shockwave.h5"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge", "vz"]
out_every    = Int(dt^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
