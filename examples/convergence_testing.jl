using PyPlot
using DelimitedFiles

derivatives = "spectral"
initial     = "Collision"

xrange = yrange = zrange = (-25., 25.)
Nx, Ny = 50, 50

m0, q0, dm, dq = 0.5, 0.0, 10.0, 3.0
dp  = (0.0, 0.0, 30.0)
x0  = (2.0, 0.0, 10.0)
sL  = (10., 10.)
sT  = (10., 10.)
mu  = dm^0.25
n   = 1

dt  = 0.025
final_time = 3.

directory    = "/home/raimon/Runs/Chihuahua/"
filename     = "test.h5"
out_file     = directory * filename
out_funcs_xz = []
out_every    = Int(dt^-1/2)

##################################

Nz_ref = 150
Nz_0   = 38
dNz    = 8

Nz_range = Nz_0:dNz:Nz_ref
func_m = []
func_q = []

for i in range(1, length(Nz_range))
    global Nz
    Nz = Nz_range[i]
    println("Now doing ", Nz)
    include("/home/raimon/dev/julia/Chihuahua/src/run.jl")
    
    push!(func_m, X[:, :, Int(Nz/2) + 1, 1])
    push!(func_q, X[:, :, Int(Nz/2) + 1, 2])
end

err_m = []
for func_i in func_m[1:end-1]
    diff_i = abs.(1 .- func_i./func_m[end])
    err_i  = maximum(diff_i)
    push!(err_m, err_i)
end

err_q = []
for func_i in func_q[1:end-1]
    diff_i = abs.(1 .- func_i./func_q[end])
    err_i  = maximum(diff_i)
    push!(err_q, err_i)
end


xx = Nz_range[1:end-1]
writedlm(directory * "convergence.csv", hcat(xx, err_m, err_q))



