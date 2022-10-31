#########################
# SOLUTION PARAMETERS
#########################

function mass()
    return X[:,:,:,1]
end

function charge()
    return X[:,:,:,2]
end

function px()
    return X[:,:,:,3]
end

function py()
    return X[:,:,:,4]
end

function pz()
    return X[:,:,:,5]
end

function vx()
    m  = X[:,:,:,1]
    p  = X[:,:,:,3]

    return (p - D(m,[1,0,0]))./m
end

function vy()
    m  = X[:,:,:,1]
    p  = X[:,:,:,4]

    return (p - D(m,[0,1,0]))./m
end

function vz()
    m  = X[:,:,:,1]
    p  = X[:,:,:,5]

    return (p - D(m,[0,0,1]))./m
end

function Px()
    m  = X[:,:,:,1]
    p  = X[:,:,:,5]

    return (p - D(m,[0,0,1]))./m
end
