using PyPlot
using LinearAlgebra
using FastGaussQuadrature

function galerkin_matrix(x)
    idx = inv.(diff(x))
    # `diff(x)` returns vector of differences,
    # i.e. `diff(x)[k] = x[k+1] - x[k]`
    return Tridiagonal(
        -idx[2:end-1],               # subdiagonal
        idx[2:end] .+ idx[1:end-1],  # diagonal
        -idx[2:end-1]                # super-diagonal
    )
end

function right_hand_side(f,x)
    dx = diff(x)
    return (dx[1:end-1] .+ dx[2:end])./2 .* f.(x[2:end-1])
end

function example()
    # Set up problem
    f = x->sin(2π*x)
    u = x->sin(2π*x)/(2π)^2
    x = LinRange(0,1,9)

    # Assemble and solve linear system
    A = galerkin_matrix(x)
    b = right_hand_side(f,x)
    c = A\b

    # Plot solutions
    clf()
    xx = LinRange(0,1,1000)
    plot(xx, u.(xx), "k-")
    # Solution is obtained by linearly interpolating the point-values stored
    # in `c` and 0 at the boundary. We use the plotting library to do the
    # interpolation for us here.
    plot(x, [0;c;0], lw=2)
end


"""
    lerp(a,b,fa,fb,x)

Evaluate the linear interpolant through `(a,fa)` and `(b,fb)` at point `x`.
"""
lerp(a,b,fa,fb, x) = ((b-x)*fa + (x-a)*fb)/(b-a)

"""
    integrate(f,a,b)

Integrate `f(x)` over the interval `[a,b]` using five-point Gauss quadrature.
"""
function integrate(f,a,b)
    qx,qw = gausslegendre(5)
    qx = (b-a)/2 .* qx .+ (b+a)/2
    return (b-a)/2 * sum(qw .* f.(qx))
end

"""
    L2_error(c,x,u)

Compute `||u - u_n||_{L^2}` with `u_n` being the finite element solution
given by mesh `x` and coefficients `c`.
"""
function L2_error(c,x,u)
    @assert length(c)+2 == length(x)
    n = length(x)-2
    e = integrate(xx->(lerp(x[1],x[2],0,c[1],xx)-u(xx))^2, x[1],x[2]) +             # first interval
        integrate(xx->(lerp(x[end-1],x[end],c[end],0,xx)-u(xx))^2, x[end-1],x[end]) # last interval
    for i = 1:n-1
        e += integrate(xx->(lerp(x[i+1],x[i+2],c[i],c[i+1],xx)-u(xx))^2, x[i+1],x[i+2])
    end
    return sqrt(e)
end

"""
    H1s_error(c,x,du)

Compute `||du/dx - du_n/dx||_{L^2}` with `u_n` being the finite element solution
given by mesh `x` and coefficients `c`.
"""
function H1s_error(c,x,du)
    @assert length(c)+2 == length(x)
    n = length(x)-2
    e = integrate(xx->(c[1]/(x[2]-x[1])-du(xx))^2, x[1],x[2]) +             # first interval
        integrate(xx->(-c[end]/(x[end]-x[end-1])-du(xx))^2, x[end-1],x[end]) # last interval
    for i = 1:n-1
        e += integrate(xx->((c[i+1]-c[i])/(x[i+2]-x[i+1])-du(xx))^2, x[i+1],x[i+2])
    end
    return sqrt(e)
end

function convergence()
    # Set up problem
    f = x->sin(2π*x)
    u = x->sin(2π*x)/(2π)^2
    du = x->cos(2π*x)/(2π)
    n = round.(Int, 10.0.^LinRange(0,3,30))

    # Compute errors
    l2 = zeros(length(n))
    h1s = zeros(length(n))
    for i = 1:length(n)
        x = LinRange(0,1,n[i]+2)
        A = galerkin_matrix(x)
        b = right_hand_side(f,x)
        c = A\b

        l2[i] = L2_error(c,x,u)
        h1s[i] = H1s_error(c,x,du)
    end

    # Plot
    clf()
    nn = (4,1e3)
    loglog(nn,inv.(nn), "k--")
    loglog(nn,1e-1.*inv.(nn).^2, "k--")
    loglog(n,l2, label="L2 error")
    loglog(n,h1s, label="H1s error")
    legend(loc="best")
end
