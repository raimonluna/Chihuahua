function initialize_outfiles()
    global funcs_array_xz, funcs_array_center, funcs_array_avg
    
    funcs_array_xz = zeros(Nx, Nz, length(out_funcs_xz), length(out_iterations_xz))
    funcs_array_center = zeros(length(out_funcs_center) + 1, Int(floor(final_time / (dt * out_every_center))) + 1)
    funcs_array_avg    = zeros(length(out_funcs_avg) + 1,    Int(floor(final_time / (dt * out_every_avg))) + 1)
end

function iteration_output()
    global funcs_array_xz, funcs_array_center, funcs_array_avg
    
    for (j, it) in enumerate(out_iterations_xz)
        if iteration == it
            for (i, func) in enumerate(out_funcs_xz)
                func_exec = @eval $(Symbol(func))
                funcs_array_xz[:, :, i, j] = func_exec()[:, Int(Ny/2) + 1, :]
            end
        end
    end
    
    if iteration % out_every_center == 0
        j = Int(round(iteration / out_every_center)) + 1
        funcs_array_center[1, j] = eff_time
        for (i, func) in enumerate(out_funcs_center)
            func_exec = @eval $(Symbol(func))
            funcs_array_center[i + 1, j] = func_exec()[Int(Nx/2) + 1, Int(Ny/2) + 1, Int(Nz/2) + 1]
        end
    end
    
    if iteration % out_every_avg == 0
    j = Int(round(iteration / out_every_avg)) + 1
        funcs_array_avg[1, j] = eff_time
        for (i, func) in enumerate(out_funcs_avg)
            func_exec = @eval $(Symbol(func))
            funcs_array_avg[i + 1, j] = Statistics.mean(func_exec())
        end
    end
    
end

function terminate_outfiles()
    global funcs_array_xz, funcs_array_center, funcs_array_avg
    
    # Remove possible existing file to avoid errors.
    try rm(out_file) catch end
    out   = h5open(out_file, "w")

    # Print the grid
    ggrid = create_group(out, "Grid")
    write_dataset(ggrid, "x", x[:,1,1])
    write_dataset(ggrid, "y", y[1,:,1])
    write_dataset(ggrid, "z", z[1,1,:])
    
    # Print the run data
    if (length(out_funcs_xz) > 0) & (length(out_iterations_xz) > 0)
        group_xz = create_group(out, "out_vars_xz")
        write_dataset(group_xz, "time", out_iterations_xz / dt)
        for (i, func) in enumerate(out_funcs_xz)
            write_dataset(group_xz, func, funcs_array_xz[:, :, i, :])
        end
    end
    
    group_center = create_group(out, "out_vars_center")
    for (i, func) in enumerate(out_funcs_center)
        write_dataset(group_center, func, funcs_array_center[[1, i + 1], :])
    end
    
    group_avg = create_group(out, "out_vars_avg")
    for (i, func) in enumerate(out_funcs_avg)
        write_dataset(group_avg, func, funcs_array_avg[[1, i + 1], :])
    end
    
    close(out)
end

