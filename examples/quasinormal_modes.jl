using PyPlot

derivatives = "spectral"
initial     = "Sound"

Nx, Ny, Nz = 2, 2, 64

m0, q0 = 1.0, 0.7
dm, dq = 0.1, (q0/m0) * 0.1

final_time = 30.0
dt  = 0.01

###############################################

directory    = "/home/raimon/Runs/Chihuahua/"
filename     = "tmp.h5"
out_file     =  directory * filename

out_funcs_xz      = []
out_iterations_xz = []

out_funcs_center = ["mass"]
out_every_center = 50

out_funcs_avg = []
out_every_avg = Inf

k_range   = collect(0.2:0.2:4.0)
all_modes = zeros(21, 61)

for i in range(1, length(k_range))
    global k, xrange, yrange, zrange
    
    k = k_range[i]
    println("Now doing k = ", k)
    xrange = yrange = zrange = (-10*pi/k, 10*pi/k)
    
    include("../src/run.jl")
    
    all_modes[1, :] = funcs_array_center[1, :]
    all_modes[i + 1, :] = funcs_array_center[2, :]
end

###############################################    

out_file = "../output/quasinormal_modes.h5"

# Remove possible existing file to avoid errors.
try rm(out_file) catch end
out   = h5open(out_file, "w")

group = create_group(out, "QNM")
write_dataset(group, "k", k_range)
write_dataset(group, "data", all_modes)

close(out)
    
    
    
    
    
