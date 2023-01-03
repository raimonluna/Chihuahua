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
initial     = "Single"

m0, q0, dm, dq = 1.0, 0.0, 20.0, 6.0
x0  = (0.0, 0.0, 10.0)
dp  = (0.0, 0.0, 60.0)
sL  = 10.
sT  = 50.
n   = 1 #For opposite charge choose n to -1

dt  = 0.1
final_time = 30.

#########################
# Output
#########################
directory    = "../output/"
filename     = "single_charged_oblate.h5"
out_file     = directory * filename

out_funcs_xz      = []
out_iterations_xz = []

out_funcs_center = ["mass", "charge", "pressure_x", "pressure_y", "pressure_z", "pressure_hydro_x", "pressure_hydro_y", "pressure_hydro_z"]
out_every_center = 1

out_funcs_avg = ["charged_entropy"]
out_every_avg = 1

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")

