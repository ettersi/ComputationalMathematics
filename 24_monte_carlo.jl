using PyPlot
using SpecialFunctions
using StaticArrays
using Statistics
using Random
using BenchmarkTools

"""
    midpoint(f,d,n)

Compute the integral of `f` over `[0,1]^d` using the midpoint rule with `n`
quadrature points in each dimension.
"""
midpoint(f,d,n) = midpoint_nested(f,n,ntuple(k->n,d))
function midpoint_nested(f,n,nn)
    q = 0.0
    x = LinRange(0,1,2n+1)[2:2:end-1]
    for i in CartesianIndices(nn)
        q += f((ik->x[ik]).(i.I))
    end
    return q/n^length(nn)
end

"""
    monte_carlo(f,d,N)

Compute the integral of `f` over `[0,1]^d` using `N` uniformly distributed
samples.
"""
monte_carlo(f,d,N) = monte_carlo(f,Val(d),N)
function monte_carlo(f,::Val{d},N) where {d}
    q = 0.0
    for i = 1:N
        q += f(@SVector rand(d))
    end
    return q/N
end

function convergence()
    d = (2,4,7,10,14)
    N = round.(Int, 10.0.^LinRange(0,6,51))
    f = x-> exp(-sum(x.^2))
    I = sqrt(π)*erf(1)/2
    ylims = [1e-8,4e0]

    close("all")
    figure(figsize=(10,4))
    # clf()

    subplot(1,2,1)
    title("Midpoint")
    for (i,d) in enumerate(d)
        n = round.(Int, N.^(1/d))
        loglog(n.^d, [abs(I^d - midpoint(f,d,n))/I^d for n in n], label="d = $d")

        NN = (1e2,N[end])
        offset = (5e-2,1e-1,2e-1,6e-1,1e0)
        loglog(NN, offset[i].*inv.(NN).^(2/d), "k--");
    end
    xlabel("N")
    ylabel("relative error")
    legend(loc="best")
    ylim(ylims)

    subplot(1,2,2)
    title("Monte Carlo")
    for d in d
        loglog(N, [abs(I^d - monte_carlo(f,d,N))/I^d for N in N])
    end
    NN = (1e2,N[end])
    loglog(NN, 6e0.*sqrt.(inv.(NN)), "k--");
    xlabel("N")
    ylabel("relative error")
    ylim(ylims)
end

function error_estimation()
    f = x-> exp(-x^2)
    I = sqrt(π)*erf(1)/2

    N = 100
    M = 1000

    F = f.(rand(M,N))
    μ = mean(F,dims=2)
    σ = std(F)

    println("Empirical root-mean-square error: ", sqrt(mean((μ.-I).^2)))
    println("Expected root-mean-square error: ", σ/sqrt(N))
    println()
    println("Empirical ratio of estimates inside one sigma confidence interval: ", sum(abs.(μ.-I) .< σ/sqrt(N))/M)
    println("Expected ratio of estimates inside one sigma confidence interval: ", erf(sqrt(0.5)))
end

function rng_benchmarks()
    prng = MersenneTwister()
    trng = RandomDevice()
    println("Pseudo-random number generator:")
    @btime rand($prng)
    println()
    println("True random number generator:")
    @btime rand($trng)
end

function transformation_sampling()
    U = rand(1_000_000)
    X = sqrt.(U)
    clf()
    hist(X; bins= 100, density = true)
end

function rejection_sampling()
    N = 1_000_000
    F = [G for G in rand(N) if rand() <= G]
    println("Acceptance ratio: ", length(F)/N)
    clf()
    hist(F; bins= 100, density = true)
end

function importance_sampling()
    # Set up problem
    # Note that g(x) is the normal probability density with mean μ and variance σ^2
    μ = 0.5
    σ = 0.01
    g = x -> 1/sqrt(2π*σ^2) * exp(-(x-μ)^2/(2*σ^2))
    h = x -> x^2 * g(x)
    f = x -> 1

    # Generate the samples
    N = 50
    F = rand(N)
    G = σ.*randn(N) .+ μ

    # Plot function and samples
    clf()
    x = LinRange(0,1,1000)
    plot(x, h.(x))
    plot(F, zero.(F), "o", ms=4)
    plot(G, zero.(G), "o", ms=4)

    # Compute errors of Monte Carlo estimates
    n = 1_000_000
    x = LinRange(0,1,2n+1)[2:2:end]
    I = sum(h.(x))/n
    println("Uniform sampling error: ", abs(sum(h.(F)/N) - I))
    println("Importance sampling error: ", abs(sum(h.(G).*f(G)./g.(G)/N) - I))
end
