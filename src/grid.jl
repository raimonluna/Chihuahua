#########################
# GENERATE GRIDS
#########################

xl  = range(xrange[1], stop=xrange[2], length=Nx+1)[1:end-1]
yl  = range(yrange[1], stop=yrange[2], length=Ny+1)[1:end-1]
zl  = range(zrange[1], stop=zrange[2], length=Nz+1)[1:end-1]

wxl = fftfreq(Nx,Nx)
wyl = fftfreq(Ny,Ny)
wzl = fftfreq(Nz,Nz)

x   = getindex.(Iterators.product(xl, yl, zl), 1) 
y   = getindex.(Iterators.product(xl, yl, zl), 2)
z   = getindex.(Iterators.product(xl, yl, zl), 3)

wx  = getindex.(Iterators.product(wxl, wyl, wzl), 1) 
wy  = getindex.(Iterators.product(wxl, wyl, wzl), 2)
wz  = getindex.(Iterators.product(wxl, wyl, wzl), 3)

#########################
# NUMERICAL METHODS
#########################

#This function calculates derivatives using the Fast Fourier Transform (FFT)
function D(f, orders)
    Dx = xrange[2] - xrange[1]
    Dy = yrange[2] - yrange[1]
    Dz = zrange[2] - zrange[1]
    
    df_fft = fft(f) .* ((im*wx).^orders[1]) .* ((im*wy).^orders[2]) .* ((im*wz).^orders[3])
    
    return real(ifft(df_fft)) * ((2*pi/Dx)^orders[1]) * ((2*pi/Dy)^orders[2]) * ((2*pi/Dz)^orders[3])
end

#Take a 4th order Runge-Kutta step 
function RK4_step!(X, dt)
    global X

    k1 = rhs(X)
    k2 = rhs(X + k1 * dt / 2)
    k3 = rhs(X + k2 * dt / 2)
    k4 = rhs(X + k3 * dt)
    
    X  = X + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
end
