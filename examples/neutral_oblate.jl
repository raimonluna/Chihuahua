#########################
# GRID
#########################

derivatives = "spectral"
xrange = yrange = (-50., 50.)
zrange = (-50.0, 50.0)
Nx, Ny, Nz = 50, 50, 150

#########################
# INITIAL DATA
#########################

#Collision_V build data with v and no initial p. 
#The variable is the same p, but it refers to v in the initial data.
initial     = "Collision"

m0, q0, dm, dq = 0.5, 0.0, 10.0, 0.0
dp  = (0.0, 0.0, 30.0)
x0  = (0.0, 0.0, 10.0)
sL  = (10., 10.)
sT  = (50., 50.)
n   = 1 #For opposite charge choose n to -1

dt  = 0.1#0.05
final_time = 30.

#########################
# Output
#########################

directory    = "/home/raimon/Downloads/"
filename     = "neutral_onlate_dx_0.h5"
out_file     = directory * filename

out_funcs_xz = ["mass", "charge"]
out_every_xz = 20

out_funcs_center = ["mass", "charge", "pressure_x", "pressure_y", "pressure_z", "pressure_hydro_x", "pressure_hydro_y", "pressure_hydro_z"]
out_every_center = 5

out_funcs_avg = ["neutral_entropy", "charged_entropy"]
out_every_avg = 5

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
