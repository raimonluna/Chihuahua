#########################
# SOLUTION PARAMETERS
#########################

function mass()
    return X[:,:,:,1]
end

function charge()
    return X[:,:,:,2]
end

#########################
# MOMENTA
#########################

function momentum_x()
    return X[:,:,:,3]
end

function momentum_y()
    return X[:,:,:,4]
end

function momentum_z()
    return X[:,:,:,5]
end

#########################
# VELOCITIES
#########################

function velocity_x()
    m  = X[:,:,:,1]
    px = X[:,:,:,3]

    return (px - D(m, [1,0,0]))./m
end

function velocity_y()
    m  = X[:,:,:,1]
    py = X[:,:,:,4]

    return (py - D(m, [0,1,0]))./m
end

function velocity_z()
    m  = X[:,:,:,1]
    pz = X[:,:,:,5]

    return (pz - D(m, [0,0,1]))./m
end

#########################
# PRESSURES
#########################

function pressure_x()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vx = velocity_x()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))
    rM    = 0.5 * (m - sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vx, [1,0,0]) .- (rP .- rM) .* D(log.(m), [2,0,0])
end

function pressure_y()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vy = velocity_y()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))
    rM    = 0.5 * (m - sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vy, [0,1,0]) .- (rP .- rM) .* D(log.(m), [0,2,0])
end

function pressure_z()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vz = velocity_z()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))
    rM    = 0.5 * (m - sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vz, [0,0,1]) .- (rP .- rM) .* D(log.(m), [0,0,2])
end

function pressure_hydro_x()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vx = velocity_x()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vx, [1,0,0])
end

function pressure_hydro_y()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vy = velocity_y()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vy, [0,1,0]) 
end

function pressure_hydro_z()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    vz = velocity_z()
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))

    return m .- rP .* D(vz, [0,0,1])
end


#########################
# ENTROPIES
#########################


function charged_entropy()
    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    
    rP    = 0.5 * (m + sqrt.(m.^2 .- 2*q.^2))

    return 4*pi*rP
end

function neutral_entropy()
    m  = X[:,:,:,1]
    vx = velocity_x()
    vy = velocity_y()
    vz = velocity_z()
    
    v2  = vx.^2 + vy.^2 + vz.^2
    dm2 = D(m, [1,0,0]).^2 + D(m, [0,1,0]).^2 + D(m, [0,0,1]).^2

    return -4*pi*(- 0.5 * m .* v2 - 0.5 * dm2 ./ m + m .* log.(m))
end










