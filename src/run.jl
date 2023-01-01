#########################
# LOAD TOOLS
#########################

using FFTW
using Printf
import Statistics
using HDF5

include("grid.jl")
include("physics.jl")
include("analysis.jl")
include("output.jl")

if derivatives == "finite differences"
    D = D4th
elseif derivatives == "spectral"
    D = DFFT
end

if initial == "Collision"
    X = Gaussian(m0, q0, dm, dq, dp, x0, sL, sT, n)
elseif initial == "Collision_V"
    X = Gaussian_V(m0, q0, dm, dq, dv, x0, sL, sT, n)
elseif initial == "Single"
    X = Single_Gaussian(m0, q0, dm, dq, dp, sL, sT)
elseif initial == "Collision_1D"
    X = Gaussian1D(m0, q0, dm, dq, dp, dz, s)
elseif initial == "Sound"
    X = Sound(k, m0, q0, dm, dq)
end

#########################
# INTEGRATION
#########################

iteration = 0
eff_time  = 0.0
toc       = Inf

initialize_outfiles()

println("----------------------------------------------------------------------")
println("Iteration   Time | et per min |      Mass density |    Charge density ")
println("                 |            | minimum   maximum | minimum   maximum ")
println("----------------------------------------------------------------------")

while eff_time <= final_time
    global iteration, eff_time, toc
    tic = time()

    @printf("%9i %6.3f | %10.3f | %7.3f %9.3f | %7.3f %9.3f \n",
    iteration, eff_time, 60*dt/toc,
    minimum(X[:,:,:,1]), maximum(X[:,:,:,1]),
    minimum(X[:,:,:,2]), maximum(X[:,:,:,2]))
    
    iteration_output()

    if eff_time < final_time
        RK4_step!(X, dt)
        iteration += 1
        eff_time  += dt
        toc = time() - tic
    end
end

terminate_outfiles()
println("Done.")
