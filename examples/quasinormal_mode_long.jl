using PyPlot

derivatives = "spectral"
initial     = "Sound"

k = 0.3
xrange = yrange = zrange = (-10*pi/k, 10*pi/k)
Nx, Ny, Nz = 2, 2, 64

m0, q0 = 1.0, 0.7
dm, dq = 0.1, (q0/m0) * 0.1

final_time = 100.0
dt  = 0.01

###############################################

directory    = "../output/"
filename     = "quasinormal_mode_long.h5"
out_file     =  directory * filename

out_funcs_xz      = []
out_iterations_xz = []

out_funcs_center = ["mass"]
out_every_center = 50

out_funcs_avg = []
out_every_avg = Inf
    
include("../src/run.jl")
    


