using LinearAlgebra
using Random
using PyPlot

function generate_data(m)
    a = rand(m)
    b = 0.2 .+ 0.6.*a .+ 0.1*randn(m)
    return a,b
end

function solve_least_squares(a,b)
    A = [ones(length(a)) a]
    Q,R = qr(A)
    return R\Matrix(Q)'*b
    # Simpler and faster: return A\b
end

function example()
    a,b = generate_data(100)
    c = solve_least_squares(a,b)

    clf()
    x = LinRange(0,1,1000)
    plot(a,b,"ko")
    plot(x, @.(c[1]+c[2]*x))
    ylim([0,1])
end
