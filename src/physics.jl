#########################
# EFFECTIVE EQUATIONS
#########################

function rhs(X)
    dXdt = zeros(Nx, Ny, Nz, 5)

    m  = X[:,:,:,1]
    q  = X[:,:,:,2]
    px = X[:,:,:,3]
    py = X[:,:,:,4]
    pz = X[:,:,:,5]

    rM    = 0.5 * (m - sqrt.(m.^2 .- 2*q.^2))

    Mxy = Myx = px.*py./m + rM .* (D(py./m, [1, 0, 0]) + D(px./m, [0, 1, 0]))
    Mxz = Mzx = px.*pz./m + rM .* (D(pz./m, [1, 0, 0]) + D(px./m, [0, 0, 1]))
    Myz = Mzy = py.*pz./m + rM .* (D(pz./m, [0, 1, 0]) + D(py./m, [0, 0, 1]))

    Mxx = px.^2 ./ m + 2 * rM.*D(px./m, [1, 0, 0])
    Myy = py.^2 ./ m + 2 * rM.*D(py./m, [0, 1, 0])
    Mzz = pz.^2 ./ m + 2 * rM.*D(pz./m, [0, 0, 1])

    dXdt[:,:,:,1] = D(m,  [2, 0, 0]) + D(m,  [0, 2, 0]) + D(m,  [0, 0, 2]) - D(px,       [1, 0, 0]) - D(py,       [0, 1, 0]) - D(pz,       [0, 0, 1])
    dXdt[:,:,:,2] = D(q,  [2, 0, 0]) + D(q,  [0, 2, 0]) + D(q,  [0, 0, 2]) - D(px.*q./m, [1, 0, 0]) - D(py.*q./m, [0, 1, 0]) - D(pz.*q./m, [0, 0, 1])

    dXdt[:,:,:,3] = D(px, [2, 0, 0]) + D(px, [0, 2, 0]) + D(px, [0, 0, 2]) - D(m, [1, 0, 0]) - D(Mxx, [1, 0, 0]) - D(Mxy, [0, 1, 0]) - D(Mxz, [0, 0, 1])
    dXdt[:,:,:,4] = D(py, [2, 0, 0]) + D(py, [0, 2, 0]) + D(py, [0, 0, 2]) - D(m, [0, 1, 0]) - D(Myx, [1, 0, 0]) - D(Myy, [0, 1, 0]) - D(Myz, [0, 0, 1])
    dXdt[:,:,:,5] = D(pz, [2, 0, 0]) + D(pz, [0, 2, 0]) + D(pz, [0, 0, 2]) - D(m, [0, 0, 1]) - D(Mzx, [1, 0, 0]) - D(Mzy, [0, 1, 0]) - D(Mzz, [0, 0, 1])

    return dXdt
end

#########################
# INITIAL DATA
#########################
#
# function GaussianSymmetric(m0, q0, dm, dq, dp, x0, s)
#
#     X = zeros(Nx, Ny, Nz, 5)
#
#     dx,  dy,  dz  = x0
#     dpx, dpy, dpz = dp
#     exp1 = exp.(- ((x .- dx).^2 + (y .- dy).^2 + (z .- dz).^2) ./ s)
#     exp2 = exp.(- ((x .+ dx).^2 + (y .+ dy).^2 + (z .+ dz).^2) ./ s)
#
#     X[:,:,:,1] = m0 .+ dm .* (exp1 + exp2)
#     X[:,:,:,2] = q0 .+ dq .* (exp1 + exp2)
#     X[:,:,:,3] = - dpx .* (exp1 - exp2)
#     X[:,:,:,4] = - dpy .* (exp1 - exp2)
#     X[:,:,:,5] = - dpz .* (exp1 - exp2)
#
#     return X
# end

function Gaussian(m0, q0, dm, dq, dp, x0, sL, sT, n)

    X = zeros(Nx, Ny, Nz, 5)

    dx,  dy,  dz  = x0
    dpx, dpy, dpz = dp
    sL1, sL2      = sL
    sT1, sT2      = sT
    exp1 = exp.(- ((x .- dx).^2 + (y .- dy).^2) ./ sT1 - (z .- dz).^2 ./ sL1)
    exp2 = exp.(- ((x .+ dx).^2 + (y .- dy).^2) ./ sT2 - (z .+ dz).^2 ./ sL2)

    X[:,:,:,1] = m0 .+ dm .* (exp1 + exp2)
    X[:,:,:,2] = q0 .+ dq .* (exp1 + n .* exp2)
    X[:,:,:,3] = - dpx .* (exp1 - exp2)
    X[:,:,:,4] = - dpy .* (exp1 - exp2)
    X[:,:,:,5] = - dpz .* (exp1 - exp2)

    return X
end

function Single_Gaussian(m0, q0, dm, dq, s)

    X    = zeros(Nx, Ny, Nz, 5)
    exp1 = exp.(- ((x).^2 + (y).^2 + (z).^2) ./ s)

    X[:,:,:,1] = m0 .+ dm .* exp1
    X[:,:,:,2] = q0 .+ dq .* exp1
    X[:,:,:,3] .= 0
    X[:,:,:,4] .= 0
    X[:,:,:,5] .= 0

    return X
end
