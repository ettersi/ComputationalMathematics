using LinearAlgebra
using FFTW
using PyPlot

function galerkin(n,f)
    x = LinRange(0,1,n+2)[2:end-1]
    A = Diagonal((π.*(1:n)).^2)
    b = FFTW.r2r(f.(x), FFTW.RODFT00) / (sqrt(2) * (n+1))
    return c = A\b
end

function evaluate(c,x)
    return sqrt(2) * sum( c[i] * sin(π*i*x) for i = 1:length(c))
end

function example()
    f = x->abs(x - 0.5) < 0.25
    c = galerkin(10,f)
    u = x->evaluate(c,x)
    x = LinRange(0,1,1000)
    clf()
    plot(x, u.(x))
end
