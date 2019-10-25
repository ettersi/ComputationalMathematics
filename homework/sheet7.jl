using PyPlot
using LinearAlgebra
using SparseArrays
using Test
using FastGaussQuadrature # install using `] add FastGaussQuadrature`

function galerkin_matrix(n)
    # TODO: your code here!
end

function right_hand_side(f,n)
    # TODO: your code here!
end

"""
   querp(a,b,fa,fm,fb,x)

Evaluate the quadratic interpolant through data points `(a,fa)`, `((a+b)/2,fm)`,
`(b,fb)` at the point `x`.
"""
function querp(a,b,fa,fm,fb,x)
    # TODO: your code here!
end

"""
   d_querp(a,b,fa,fm,fb,x)

Evaluate the derivative of the quadratic interpolant through data points
`(a,fa)`, `((a+b)/2,fm)`, `(b,fb)` at the point `x`.
"""
function d_querp(a,b,fa,fm,fb,x)
    # TODO: your code here!
end

function test_querp()
    @testset "querp" begin
        @test querp(0,2, 0,1,2, 0.1) ≈ 0.1
        @test querp(0,2, 0,1,4, 0.1) ≈ 0.01
    end
end

function test_d_querp()
    @testset "d_querp" begin
        @test d_querp(0,2, 0,1,2, 0.1) ≈ 1
        @test d_querp(0,2, 0,1,4, 0.1) ≈ 0.2
    end
end

function evaluate_solution_on_fine_mesh(c,nn)
    n = (length(c)-1)÷2
    x = LinRange(0,1,(n+1)*nn+1)
    u = zeros((n+1)*nn+1)
    @. u[1:nn] = querp(x[1],x[nn+1], 0,c[1],c[2], x[1:nn])
    for i = 1:n-1
        @. u[i*nn + (1:nn)] = querp(
            x[i*nn+1], x[(i+1)*nn+1],
            c[2i], c[2i+1], c[2i+2],
            x[i*nn+1:(i+1)*nn]
        )
    end
    @. u[end-nn:end] = querp(x[end-nn],x[end], c[end-1],c[end],0, x[end-nn:end])
    return x,u
end

function plot_solution()
    # Set up problem
    f = x->sin(2π*x)
    u = x->sin(2π*x)/(2π)^2
    n = 3

    # Assemble and solve linear system
    A = galerkin_matrix(n)
    b = right_hand_side(f,n)
    c = A\b

    # Plot reference solution
    clf()
    x = LinRange(0,1,1000)
    plot(x, u.(x), "k-")

    # Plot Galerkin solution
    x,u = evaluate_solution_on_fine_mesh(c,20)
    plot(x,u)
end


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
function L2_error(c,u)
    n = (length(c)-1)÷2
    x = LinRange(0,1,n+2)
    e = integrate(xx->(querp(x[1],x[2],0,c[1],c[2],xx)-u(xx))^2, x[1],x[2]) +           # first interval
        integrate(xx->(querp(x[end-1],x[end],c[end-1],c[end],0,xx)-u(xx))^2, x[end-1],x[end]) # last interval
    for i = 1:n-1
        e += integrate(xx->(querp(x[i+1],x[i+2],c[2i],c[2i+1],c[2i+2],xx)-u(xx))^2, x[i+1],x[i+2])
    end
    return sqrt(e)
end

"""
    H1s_error(c,x,du)

Compute `||du/dx - du_n/dx||_{L^2}` with `u_n` being the finite element solution
given by mesh `x` and coefficients `c`.
"""
function H1s_error(c,du)
    n = (length(c)-1)÷2
    x = LinRange(0,1,n+2)
    e = integrate(xx->(d_querp(x[1],x[2],0,c[1],c[2],xx)-du(xx))^2, x[1],x[2]) +           # first interval
        integrate(xx->(d_querp(x[end-1],x[end],c[end-1],c[end],0,xx)-du(xx))^2, x[end-1],x[end]) # last interval
    for i = 1:n-1
        e += integrate(xx->(d_querp(x[i+1],x[i+2],c[2i],c[2i+1],c[2i+2],xx)-du(xx))^2, x[i+1],x[i+2])
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
        A = galerkin_matrix(n[i])
        b = right_hand_side(f,n[i])
        c = A\b

        l2[i] = L2_error(c,u)
        h1s[i] = H1s_error(c,du)
    end

    # Plot
    clf()
    nn = (4,1e3)
    loglog(nn,inv.(nn).^2, "k--")
    loglog(nn,1e-1.*inv.(nn).^3, "k--")
    loglog(n,l2, label="L2 error")
    loglog(n,h1s, label="H1s error")
    legend(loc="best")
end
