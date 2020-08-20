using PyPlot
using StatsBase
using Random

function rounding_error()
    Random.seed!(42)  # Make `randn()` reproducible

    n = 100_000_000
    a = randn(n)
    b = randn(n)

    c64 = @. Float64(a) + Float64(b)
    c32 = @. Float32(a) + Float32(b)
    relative_error = @. abs(c64-c32)/abs(c64)

    log10_x = LinRange(-14,0,100)
    x = 10.0.^log10_x
    m = (x[2:end] .+ x[1:end-1]) ./ 2
    h = fit(Histogram, log10.(relative_error), log10_x).weights ./ n

    clf()
    loglog(m, h)
    loglog(eps(Float32).*[1,1], [1e-7,1e-1], "k--", label="eps(Float32)")
    xlabel("rounding error")
    ylabel("frequency")
    legend(loc="best", frameon=false)
    tight_layout()
    display(gcf())

    println()
    println("P(relative_error < 1e-4) = ", sum(relative_error .< 1e-4)/n)
    println("P(relative_error > 1e-1) = ", sum(relative_error .> 1e-1)/n)
end


function condition_number()
    f = sin; df = cos
    x = 1.0; x̃ = x + 1e-8

    exact_relative_error = abs( ( f(x̃) - f(x) ) / f(x) )

    κ = abs( df(x)*x / f(x) )
    estimated_relative_error = κ * abs( (x̃ - x) / x )

    println("    Exact relative error: ", exact_relative_error)
    println("Estimated relative error: ", estimated_relative_error)
end


function stability()
    x = 10.0.^LinRange(-16,0,100)
    clf()
    for (label,f) = (
        ("unstable", x -> sqrt(1 + x) - 1),
        (  "stable", x -> x / (sqrt(1+x) + 1))
    )
        relative_error = @. abs(f(x) .- f(big(x))) / abs(f(big(x)))
        loglog(x, relative_error; label=label)
    end
    xlabel(L"x")
    ylabel("Relative error")
    legend(loc="best", frameon=false)
    display(gcf())
end