# Linear algebra packages
using LinearAlgebra
using SparseArrays
using Random
using IterativeSolvers
using Preconditioners
# If any of the above lines throw an error, it is most likely because  you have
# not installed the corresponding package. To install it, type
#    ] add PackageName
# on the REPL


# Plotting packages
using PyPlot
using LaTeXStrings
ioff() # Disable interactive mode. It messes with small figsizes
rcdefaults()
rc("text", usetex=true)
rc("text.latex", preamble="\\usepackage{amsmath},\\usepackage{amssymb},\\usepackage{sfmath}")
rc("font", size=8)


laplacian_1d(n) = Tridiagonal(
    fill(Float64((n+1)^2),n-1), # subdiagonal
    fill(Float64(-2*(n+1)^2),n),   # diagonal
    fill(Float64((n+1)^2),n-1)  # superdiagonal
)

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,(n,n))
    return kron(Δ,Id) + kron(Id,Δ)
end


function well_conditioned_matrix()
    κ = 50
    A = Diagonal(LinRange(1,κ,1000))
    b = rand(1000)

    n = 50
    fig = figure(figsize=(2.4,1.8))
    for (name,method) in (("MinRes",minres), ("CG",cg))
        _,ch = method(A,b, log=true, maxiter=n)
        semilogy(0:length(ch[:resnorm])-1, ch[:resnorm]./ch[:resnorm,1], label=name)
    end
    nn = [5,n]
    ρ = (sqrt(κ) - 1) / (sqrt(κ) + 1)
    semilogy(nn, 4*ρ.^nn, "k--", label=L"\mathcal{O}\bigl(\rho^{n}\bigr)")
    legend(loc="best", frameon=false)
    xlabel(L"n")
    ylabel(L"\|r_n\|_2 / \|r_0\|_2")
    yticks(10.0.^(-6:2:1))
    tight_layout()
    savefig("pics/well_conditioned_matrix.pdf")
    close(fig)
end


function one_dimensional_laplacian()
    A = -laplacian_1d(1000)
    b = rand(1000)

    n = 1000
    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    for (name,method) in (("MinRes",minres), ("CG",cg))
        _,ch = method(A,b, log=true, maxiter=n)
        semilogy(0:length(ch[:resnorm])-1, ch[:resnorm]./ch[:resnorm,1], label=name)
    end
    ylim([1e-8,4e1])
    legend(loc="best", frameon=false)
    xlabel(L"n")
    ylabel(L"\|r_n\|_2 / \|r_0\|_2")
    yticks(10.0.^(-8:2:2))

    subplot(1,2,2)
    n = 10:10:1000
    niters = [cg(laplacian_1d(n),rand(n), log=true, tol=1e-8)[2].iters for n in n]
    loglog(n, niters, "C2")
    loglog(n, 2*n, "k--", label=L"\mathcal{O}\bigl(N\bigr)")
    xlabel(L"N")
    ylabel(L"n")
    legend(loc="best", frameon=false)

    tight_layout()
    savefig("pics/one_dimensional_laplacian.pdf")
    close(fig)
end

function two_dimensional_laplacian()
    A = -laplacian_2d(100)
    b = rand(100^2)

    n = 250
    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    for (name,method) in (("MinRes",minres), ("CG",cg))
        _,ch = method(A,b, log=true, maxiter=n)
        semilogy(0:length(ch[:resnorm])-1, ch[:resnorm]./ch[:resnorm,1], label=name)
    end
    nn = [0,n]
    κ = 4*101^2/π^2
    semilogy(nn, 8*((sqrt(κ) - 1) / (sqrt(κ) + 1)).^nn, "k--")
    xlabel(L"n")
    ylabel(L"\|r_n\|_2 / \|r_0\|_2")
    yticks(10.0.^(-6:2:0))
    legend(loc="lower left", frameon=false)

    subplot(1,2,2)
    n = 10:10:100
    niters = [cg(laplacian_2d(n),rand(n^2), log=true, tol=1e-8)[2].iters for n in n]
    loglog(n.^2, niters, "C2")
    loglog(n.^2, 4n, "k--", label=L"\mathcal{O}\bigl(N^{1/2}\bigr)")
    legend(loc="best", frameon=false)
    xlabel(L"N")
    ylabel(L"n")

    tight_layout()
    savefig("pics/two_dimensional_laplacian.pdf")
    close(fig)
end


function two_dimensional_laplacian_with_ilu()
    A = -laplacian_2d(100)
    b = rand(100^2)

    n = 80
    fig = figure(figsize=(4.2,1.6))
    subplot(1,2,1)
    _,ch = cg(A,b, log=true, maxiter=n, Pl = CholeskyPreconditioner(A,0))
    semilogy(0:length(ch[:resnorm])-1, ch[:resnorm]./ch[:resnorm,1], "C1", label="CG")

    nn = [0,n]
    κ = 4*101^2/π^2
    semilogy(nn, 8*((sqrt(κ) - 1) / (sqrt(κ) + 1)).^(2 .* nn), "k--", label=L"\mathcal{O}\bigl(\rho^{n}\bigr)")
    xlabel(L"n")
    ylabel(L"\|r_n\|_2 / \|r_0\|_2")

    Random.seed!(42)
    subplot(1,2,2)
    n = 10:10:100
    niters = [begin
        A = -laplacian_2d(n)
        b = rand(n^2)
        _,ch = cg(A,b, log=true, tol=1e-8, Pl = CholeskyPreconditioner(A,0))
        ch.iters
    end for n in n]
    loglog(n.^2, niters, "C2")
    loglog(n.^2, 2n, "k--", label=L"\mathcal{O}\bigl(N^{1/2}\bigr)")
    legend(loc="best", frameon=false)
    xlabel(L"N")
    ylabel(L"n")

    tight_layout()
    savefig("pics/two_dimensional_laplacian_with_ilu.pdf")
    close(fig)
end
