using LinearAlgebra
using PyPlot

poisson(N) = Tridiagonal(
    fill(-(N+1)^2,N-1),
    fill(2*(N+1)^2,N),
    fill(-(N+1)^2,N-1)
)

##########################
# Jacobi and Gauss-Seidel

function jacobi_step(x0,b)
    N = length(x0)
    x1 = zeros(N)

    x1[1] = 0.5*(b[1]/(N+1)^2 + x0[2])
    for i = 2:N-1
        x1[i] = 0.5*(b[i]/(N+1)^2 + x0[i-1] + x0[i+1])
    end
    x1[N] = 0.5*(b[N]/(N+1)^2 + x0[N-1])

    return x1
end

function gauss_seidel_step(x,b)
    N = length(x)

    x[1] = 0.5*(b[1]/(N+1)^2 + x[2])
    for i = 2:N-1
        x[i] = 0.5*(b[i]/(N+1)^2 + x[i-1] + x[i+1])
    end
    x[N] = 0.5*(b[N]/(N+1)^2 + x[N-1])

    return x
end

function convergence_history(b, step, kmax)
    rhist = zeros(kmax)
    x = copy(b)
    for k = 1:kmax
        x = step(x,b)
        rhist[k] = norm(b - poisson(length(x))*x)
    end
    return rhist
end

function plot_convergence()
    N = 127
    kmax = 10_000
    b = ones(N)

    clf()
    for (name,step) in (
            ("Jacobi", jacobi_step),
            ("Gauss-Seidel",gauss_seidel_step),
            # ("Multigrid",multigrid_step),
        )
        rhist = convergence_history(b, step, kmax)
        semilogy(1:kmax, rhist./rhist[1], label=name)
    end

    λmax = cos(π / (N+1))
    semilogy(λmax.^(1:kmax), "k--")

    legend(loc="best", frameon=false)
end

function speed_of_propagation()
    N = 6
    b = zeros(N); b[1] = (N+1)^2
    x = 1 .- (1:N)./(N+1)
    clf()
    plot(1:N, x,"k", label="exact")
    for (i,k) in enumerate((1,2,3))
        xk = b/(N+1)^2
        for kk = 1:k
            xk = jacobi_step(xk,b)
        end
        plot(1:N, xk, label="k = $k")
    end
    legend(loc="best", frameon=false)
end



############
# Multigrid

function relaxed_jacobi_step(x0,b)
    N = length(x0)
    x1 = zeros(N)

    x1[1] = x0[1]/3 + 2/3 * 0.5*(b[1]/(N+1)^2 + x0[2])
    for i = 2:N-1
        x1[i] = x0[i]/3 + 2/3 * 0.5*(b[i]/(N+1)^2 + x0[i-1] + x0[i+1])
    end
    x1[N] = x0[N]/3 + 2/3 * 0.5*(b[N]/(N+1)^2 + x0[N-1])

    return x1
end

approximate(x) = @.(0.25 * x[1:2:end-2] + 0.5 * x[2:2:end-1] + 0.25 * x[3:2:end])

function interpolate(xx)
    x = zeros(2*length(xx) + 1)
    x[1] = 0.5*xx[1]
    @. x[2:2:end] = xx
    @. x[3:2:end-2] = 0.5 * (xx[1:end-1] + xx[2:end])
    x[end] = 0.5*xx[end]
    return x
end

function twogrid_step(x0,b)
    N = length(x0)
    @assert isodd(N)
    A = poisson(N)
    NN = N÷2 # Integer division. Rounds down if N is odd.
    AA = poisson(NN)

    x̃0 = relaxed_jacobi_step(x0,b)
    r = b - A*x̃0
    rr = approximate(r)
    Δxx = AA \ rr
    Δx = interpolate(Δxx)
    x̃1 = x̃0 + Δx
    x1 = relaxed_jacobi_step(x̃1,b)

    return x1
end

function multigrid_step(x0,b)
    N = length(x0)
    @assert isodd(N)
    A = poisson(N)
    NN = N÷2 # Integer division. Rounds down if N is odd.

    if N < 8
        # Problem is small enough that we can solve it using LU
        return A \ b
    end

    x̃0 = relaxed_jacobi_step(x0,b)
    r = b - A*x̃0
    rr = approximate(r)
    Δxx = multigrid_step(zeros(NN), rr)
    Δx = interpolate(Δxx)
    x̃1 = x̃0 + Δx
    x1 = relaxed_jacobi_step(x̃1,b)

    return x1
end

function plot_multigrid_convergence()
    kmax = 15

    clf()
    for (line,step) in (
        ("-.",twogrid_step),
        ("-",multigrid_step),
    )
        for (i,N) = enumerate((127,255))
            b = ones(N)
            rhist = convergence_history(b, step, kmax)
            semilogy(rhist./rhist[1], line*"C$(i-1)", label="N = $N")
        end
    end
    semilogy(1e-1*(1/9).^(0:kmax-1), "k--")
    legend(loc="best", frameon=false)
end
