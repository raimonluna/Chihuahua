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

m0, q0, dm, dq = 1.0, 0.0, 20.0, 0.0
dp  = (0.0, 0.0, 60.0)
sL  = (10., 10.)
sT  = (50., 50.)
n   = 1 #For opposite charge choose n to -1

dt  = 0.1
final_time = 30.

#########################
# Output
#########################
directory    = "../output/"

out_funcs_center = ["mass", "pressure_x", "pressure_y", "pressure_z", "pressure_hydro_x", "pressure_hydro_y", "pressure_hydro_z"]
out_every_center = 1

out_funcs_avg = ["neutral_entropy"]
out_every_avg = 1

#########################
# LAUNCH RUN
#########################

for dx in [0, 2, 5, 8, 10]
    global x0, out_file, out_funcs_xz, out_iterations_xz
    
    if dx == 2
        out_funcs_xz      = ["mass"]
        out_iterations_xz = [15]
    else
        out_funcs_xz      = []
        out_iterations_xz = []
    end

    filename     = "neutral_oblate_dx_" * string(dx) * ".h5"
    out_file     = directory * filename
    x0  = (dx, 0.0, 10.0)
    
    include("../src/run.jl")
end

