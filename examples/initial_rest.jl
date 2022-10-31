#########################
# PARAMETERS
#########################
derivatives = "spectral"
initial     = "Collision"

xrange = yrange = (-30., 30.)
zrange = (-30.0, 30.0)
Nx, Ny, Nz = 48, 48, 48

m0, q0, dm, dq = 0.5, 0.1, 50.0, 15.0
dp  = (0.0, 0.0, 0.0)
x0  = (0.0, 0.0, 10.0)
sL  = (4., 4.)
sT  = (4., 4.)
mu  = dm^0.25

dt = 0.1 / mu
final_time = 30 / mu + dt

directory    = "/Users/mikel.sanchez/Documents/Chihuahua/data/initial_rest/charged_brane/"
filename     = "S_dz_10_m0_0.5_dm_50_q0_0.1_dq_15_L_60_Nx_Ny_48_Nz_48_t_30_dt_0.1.m"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
