#########################
# PARAMETERS
#########################
derivatives = "spectral"
initial     = "Single"

xrange = yrange = (-30., 30.)
zrange = (-30., 30.)
Nx, Ny, Nz = 48, 48, 48

m0, q0, dm, dq = 0.5, 0.1, 50.0, 15.0
s   = 4.
mu  = dm^0.25

dt = 0.1 / mu
final_time = 30 / mu + dt

directory    = "/Users/mikel.sanchez/Documents/Chihuahua/data/single/charged_brane/"
filename     = "S_L_60_N_48_dm_50_dq_15_m0_0.5_q0_0.1_t_30_dt_0.1.m"
out_file     = directory * filename
out_funcs_xz = ["mass", "charge"]
out_every    = Int((dt * mu)^-1)

#########################
# LAUNCH RUN
#########################

include("../src/run.jl")
