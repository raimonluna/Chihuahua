#########################
# LOAD TOOLS
#########################

using FFTW
using Printf

include("grid.jl")
include("physics.jl")
include("analysis.jl")
include("utils.jl")

#########################
# INTEGRATION
#########################

iteration = 0
eff_time  = 0.0
toc       = Inf

out = open(out_file, "w")
X = GaussianSymmetric(m0, q0, dm, dq, dp, x0, s)

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

    if iteration % out_every == 0
        for func in out_funcs_xz
            func_exec = @eval $(Symbol(func))
            write(out, ExportToMathematicaInterp(func_exec()[:, Int(Ny/2) + 1, :, 1], func * "It" * string(iteration)))
        end
    end
    
    RK4_step!(X, dt)
    iteration += 1
    eff_time  += dt
    
    toc = time() - tic
end

close(out)
println("Done.")

