using LinearAlgebra
using SparseArrays
using PyPlot

laplacian(n) = n^2 * sparse(Tridiagonal(
    fill(-1.0,n-1),
    fill( 2.0,n),
    fill(-1.0,n-1)
))

function first_order_matrix(n)
    A = laplacian(n)
    A[n,n-1:n] .= n.*[-1,1]
    return A
end

function second_order_matrix(n)
    A = laplacian(n)
    A[n,n-2:n] .= n.*[0.5,-2.0,1.5]
    return A
end

function errors(mat,n)
    u = x->sin(π/2*x)
    f = x->π^2/4 * sin(π/2*x)
    return [begin
        x = LinRange(0,1,n+1)[2:end]
        A = mat(n)
        b = [f.(x[1:end-1]); 0]
        ũ = A\b
        norm(ũ .- u.(x), 2)/sqrt(n)
    end for n in n]
end

function convergence()
    n = 2 .^ (2:12)
    first_order_errors = errors(first_order_matrix,n)
    second_order_errors = errors(second_order_matrix,n)

    fig = figure(figsize=(4,3))
    loglog(n,first_order_errors, label="error u_n")
    loglog(n,second_order_errors, label="error û_n")
    nn = 2 .^ (6:12)
    loglog(nn, 1.5*nn.^-1, "k--", label="O(n^-1)")
    loglog(nn, 0.25*nn.^-2, "k-.", label="O(n^-2)")
    xlabel("n")
    legend(loc="best")
    savefig("sheet3.png")
    close(fig)
end
