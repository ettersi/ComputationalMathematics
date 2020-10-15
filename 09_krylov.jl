using PyPlot
using LinearAlgebra
using FFTW
using BernsteinEllipses

chebyshev_points(n) = cos.(LinRange(0,π,2n+1)[2:2:end-1])

"""
    chebyshev_coefficients(f,n) -> c

Compute coefficients `c[k]` such that

    p(x) = sum( c[k] * T_k(x) for k = 0:n-1 )

is the `n`-point Chebyshev interpolant to `f(x)`.
`T_k(x)` denotes the `k`th Chebyshev polynomial.
"""
function chebyshev_coefficients(f,n)
    # WARNING:
    # This function uses some advanced techniques to compute the Chebyshev
    # coefficients efficiently and  conveniently. You are not expected to
    # understand how this function works.

    x = chebyshev_points(n)
    c = FFTW.r2r(f.(x), FFTW.REDFT10) ./ n
    c[1] /= 2
    return c
end

"""
    chebyshev_sum(c,x)

Evaluate `sum( c[k] * T_k(x) for k = 0:n-1 )` where `T_k(x)` denotes the
`k`th Chebyshev polynomial.
"""
function chebyshev_sum(c,x)
    T = (1.0,x)
    p = c[1]*T[1] + c[2]*T[2]
    for i = 3:length(c)
        T = T[2], 2*x*T[2] - T[1]
        p += c[i]*T[2]
    end
    return p
end

"""
    chebyshev_sum(c,A,b)

Evaluate `sum( c[k] * T_k(A) * b for k = 0:n-1 )` where `T_k(x)` denotes the
`k`th Chebyshev polynomial.
"""
function chebyshev_sum(c,A,b)
    T = (b,A*b)
    p = c[1]*T[1] + c[2]*T[2]
    for i = 3:length(c)
        T = T[2], 2*A*T[2] - T[1]
        p += c[i]*T[2]
    end
    return p
end

"""
    solve(A,b,m,c,d)

Evaluate `p(A)*b` where `p` is the `m`-point Chebyshev interpolant to
`x->inv(x)`on the interval `[c,d]`.
"""
function solve(A,b,m,c,d)
    # Linear map ψ : [-1,1] -> [c,d] and its inverse
    # Multiplying constant terms by `I` makes sure the function works for both
    # scalar and matrix arguments
    ψ = x̂ -> (d+c)/2*I + (d-c)/2*x̂
    ψinv = x -> (x - (d+c)/2*I) * 2/(d-c)

    coeffs = chebyshev_coefficients(inv ∘ ψ, m)
    return chebyshev_sum(coeffs, ψinv(A), b)
end

function example()
    n = 100     # size of the linear system
    m = 10      # degree of polynomial approximation
    c,d = 1,2   # Interval of eigenvalues

    # Assemble and solve linear system
    λ = LinRange(Sa,Sb,n)
    A = Diagonal(λ)
    b = ones(n)
    pAb = solve(A,b,m,c,d)

    # Assemble polynomial for comparison
    ψ = x -> (d+c)/2 + (d-c)/2*x
    ψinv = x -> (x - (d+c)/2) * 2/(d-c)
    coeffs = chebyshev_coefficients(inv ∘ ψ, m)
    p = x -> chebyshev_sum(coeffs,ψinv(x))

    # Print errors
    bound = norm(inv.(λ) .- p.(λ), Inf) * norm(b)
    error = norm(A\b .- pAb,2)
    println(" Error bound: ", bound)
    println("Actual error: ", error)
end

function convergence()
    m = 2:20
    n = 100
    κ = 10

    errors = [begin
        λ = LinRange(1,κ,n)
        A = Diagonal(λ)
        b = ones(n)
        pAb = solve(A,b,m,1,κ)
        norm(A\b .- pAb,2)
    end for m in m]

    clf()
    mm = (10,20)
    semilogy(m,errors)
    semilogy(mm, 1e1.*((sqrt(κ)-1)/(sqrt(κ)+1)).^mm, "k--")
    xlabel(L"m")
    ylabel("error")
    display(gcf())
end

function polynomial()
    c,d = 1,2        # Interval of eigenvalues
    m = 34,38,50     # degrees of polynomial approximation
    xlims = (0.5,2.5)
    ylims = (0,2)

    clf()
    x = LinRange(xlims[1], xlims[2], 1000)
    plot(x, inv.(x), "k-", label=L"1/x")
    for m in m
        # Assemble polynomial approximation to `inv`
        ψ = x -> (d+c)/2 + (d-c)/2*x
        ψinv = x -> (x - (d+c)/2) * 2/(d-c)
        coeffs = chebyshev_coefficients(inv ∘ ψ, m)
        p = x -> chebyshev_sum(coeffs,ψinv(x))

        plot(x, p.(x), label="m = $m")
    end
    fill_between([xlims[1],c],ylims[1].*[1,1],ylims[2].*[1,1], color="0.9")
    fill_between([d,xlims[2]],ylims[1].*[1,1],ylims[2].*[1,1], color="0.9")
    xlim(xlims)
    ylim(ylims)
    xlabel(L"x")
    legend(frameon=false, loc="upper center")
    display(gcf())
end

