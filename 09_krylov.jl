using PyPlot
using LinearAlgebra
using FFTW
using IterativeSolvers
using Random
using SparseArrays

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
    λ = LinRange(c,d,n)
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

"""
    krylov_vectors(A,b,m) -> V

Assemble a matrix `V` such that `V[:,k] = A^(k-1) * b`.
"""
function krylov_vectors(A,b,m)
    n = length(b)
    @assert size(A) == (n,n)

    V = zeros(n,m)
    for l = 1:m
        V[:,l] = b
        b = A*b
    end
    return V
end

function gmres_unstable(A,b,m)
    V = krylov_vectors(A,b,m)
    y = ((A*V) \ b)
    return V * y
end


function gmres_convergence()
    κ = 3
    n = 100
    m = 1:25
    A = Diagonal(LinRange(1,κ,n))
    x = ones(n)
    errors = [begin
        if (unstable = true)
            x̃ = gmres_unstable(A,A*x,m)
        else
            x̃ = gmres(A,A*x, tol = 0, maxiter = m, restart = m)
        end
        norm(x̃ - x,2)
    end for m in m]

    clf()
    semilogy(m, errors, "o-")
    mm = 1:10
    semilogy(mm, ((sqrt(κ)-1)/(sqrt(κ)+1)).^mm, "k--",
        label=L"O\left(\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^m\right)")
    xlabel(L"\mathrm{degree}(p)")
    ylabel(L"\|p(A) \, b - A^{-1} b\|_2")
    legend(frameon=false)
    display(gcf())
end

function restarted_gmres_good()
    n = 100
    A = Diagonal(LinRange(1,10,n))
    b = ones(n)

    clf()
    for (i,k) = enumerate((2,5,10,100))
        _,log = gmres(A,b; log=true, restart = k)
        semilogy(
            1:log.iters,
            log[:resnorm],
            "C$(i-1)-"
        )
        semilogy(
            1:k:log.iters,
            log[:resnorm][1:k:end],
            "C$(i-1)o",
            ms = 4
        )
        semilogy(
            [NaN], [NaN],
            "C$(i-1)-o",
            label="restart = $k",
            ms = 4
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end

function restarted_gmres_bad()
    n = 100
    A = Diagonal([0.1; LinRange(1,10,n-1)])
    b = ones(n)

    clf()
    for (i,k) = enumerate((2,5,10,100))
        _,log = gmres(A,b; log=true, restart = k)
        semilogy(
            1:log.iters,
            log[:resnorm],
            "C$(i-1)-"
        )
        semilogy(
            1:k:log.iters,
            log[:resnorm][1:k:end],
            "C$(i-1)o",
            ms = 4
        )
        semilogy(
            [NaN], [NaN],
            "C$(i-1)-o",
            label="restart = $k",
            ms = 4
        )
    end
    xlabel(L"k")
    ylabel(L"\|Ax_k - b\|_2")
    legend(loc="best")
    display(gcf())
end



function gmres_vs_minres()
    n = 200
    Random.seed!(42)
    A = rand(n,n)
    A = A+A' + 12*I
    b = rand(n)

    clf()
    for (label, log) = (
        ("GMRES", gmres(A,b, log=true, restart=length(b))[2]),
        ("MinRes", minres(A,b, log=true)[2]),
        # ("GMRES(11)", gmres(A,b, log=true, restart=11)[2]),
    )
        semilogy(log[:resnorm], "-o", ms=2, label=label)
    end
    xlabel(L"m")
    ylabel(L"\|A \, p(A) \, b - b\|_2")
    legend(frameon=false)
    display(gcf())
end



function laplacian_1d(n)
    return (n+1)^2 * Tridiagonal(
        fill( 1.0,n-1), # subdiagonal
        fill(-2.0,n),   # diagonal
        fill( 1.0,n-1)  # superdiagonal
    )
end

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Id,Δ) + kron(Δ,Id)
end

function cg_poisson_1d()
    clf()
    for n = (500,1000,1500,2000)
        Random.seed!(42)
        A = -laplacian_1d(n)
        b = rand(n)
        r = cg(A,b, log = true, tol = eps())[2][:resnorm]
        semilogy(0:length(r)-1, r, label=latexstring("n = $n"))
    end
    xlabel(L"m")
    ylabel(L"\|A \, p(A) \, b - b\|_2")
    legend(frameon=false)
    display(gcf())
end

function cg_poisson_2d()
    clf()
    for n = (50,100,150,200)
        Random.seed!(42)
        κ = 4*(n+1)^2/π^2
        A = -laplacian_2d(n)
        b = rand(n^2)
        r = cg(A,b, log = true, tol = eps())[2][:resnorm]
        semilogy(0:length(r)-1, r, label=latexstring("n = $n"))
    end
    xlabel(L"m")
    ylabel(L"\|A \, p(A) \, b - b\|_2")
    legend(frameon=false)
    display(gcf())
end