#########################
# PARAMETERS
#########################
derivatives = "finite differences"

xrange = yrange = zrange = (-40.0, 40.0)
Nx, Ny, Nz = 150, 150, 150

m0, q0, dm, dq = 0.5, 0.0, 10.0, 6.0
dp  = (0.0, 0.0, 30.0)
x0  = (5.0, 0.0, 10.0)
sL  = (10., 10.)
sT  = (50., 50.)
mu  = dm^0.25

dt = 0.05 / mu
final_time = 20 / mu + dt

directory    = "/home/mikel/Dropbox/Projects/Large D shockwaves/Data_Chihuahua/Comparison_Mathematica/"
filename     = "FD_dx_5_L_80_Nx_Ny_150_Nz_150_t_20_dt_0.05.m"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
