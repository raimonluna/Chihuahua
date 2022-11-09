#########################
# LOAD TOOLS
#########################

using FFTW
using Printf
using HDF5

include("grid.jl")
include("physics.jl")
include("analysis.jl")
include("utils.jl")

if derivatives == "finite differences"
    D = D4th
elseif derivatives == "spectral"
    D = DFFT
end

#########################
# INTEGRATION
#########################

iteration = 0
eff_time  = 0.0
toc       = Inf

# out = open(out_file, "w")
# Remove possible existing file to avoid errors.
try rm(out_file) catch end
out   = h5open(out_file, "w")
# Print the grid
ggrid = g_create(out, "Grid")
write_dataset(ggrid, "x", x[:,1,1])
write_dataset(ggrid, "y", y[1,:,1])
write_dataset(ggrid, "z", z[1,1,:])
write_dataset(ggrid, "t", Array(0.:dt*out_every:final_time))
write_dataset(ggrid, "dIt", Int(floor(out_every)))
write_dataset(ggrid, "Itmax", Int(floor(final_time/dt-1)))
write_dataset(ggrid, "mu", mu)
write_dataset(ggrid, "m0", m0)
write_dataset(ggrid, "q0", q0)

if initial == "Collision"
    X = Gaussian(m0, q0, dm, dq, dp, x0, sL, sT, n)
elseif initial == "Collision_V"
    X = Gaussian_V(m0, q0, dm, dq, dv, x0, sL, sT, n)
elseif initial == "Single"
    X = Single_Gaussian(m0, q0, dm, dq, s)
elseif initial == "Collision_1D"
    X = Gaussian1D(m0, q0, dm, dq, dp, dz, s)
end


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

    # if iteration % out_every == 0
    #     for func in out_funcs_xz
    #         func_exec = @eval $(Symbol(func))
    #         write(out, ExportToMathematicaInterp(func_exec()[:, Int(Ny/2) + 1, :, 1], func * "It" * string(iteration)))
    #     end
    # end
    #HDF5 output
    if iteration % out_every == 0
        group = g_create(out, "It = " * string(iteration, pad = 3))
        for func in out_funcs_xz
            func_exec = @eval $(Symbol(func))
            write_dataset(group, func, func_exec())
        end
    end

    RK4_step!(X, dt)
    iteration += 1
    eff_time  += dt

    toc = time() - tic
end

close(out)
println("Done.")
