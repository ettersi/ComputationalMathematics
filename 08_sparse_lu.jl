using PyPlot
using LinearAlgebra
using SparseArrays
using AMD

using Printf

function wilkinson(n)
    A = I - tril(ones(n,n),-1)
    A[:,end] .= 1
    return A
end

function wilkinson()
    n = 100
    A = wilkinson(n)
    x = ones(n)
    x̃ = A\(A*x)

    println()
    println("Condition number:")
    # You have shown in Assignment 1, Question 4 that κ(A ↦ A⁻¹b) = κ(A)
    println("    ", cond(A))
    println()
    println("Expected relative error for a backward-stable algorithm:")
    println("    ", cond(A)*eps())
    println()
    println("Empirical relative error for LU factorisation:")
    println("    ", norm(x-x̃)/norm(x))
    println()
end


function time_lse(n)
    A = rand(n,n)
    b = rand(n)
    @elapsed A\b
end

function time_lse()
    time_lse(1) # Dummy execution
    # Julia compiles functions the first time they are run.
    # To avoid compilation interfering with the benchmark,
    # we execute the function once and discard its result.

    println()
    for n = 2 .^ (6:13)
        @printf("n = %4u: %.3f seconds\n", n, time_lse(n))
    end
end



using Random

function fill_in_example()
    Random.seed!(42) # Change the argument to obtain different random matrices

    A = Matrix(5I + sprand(5,5,0.5))
    println("A:")
    display(A .!= 0); println()
    L,U = lu(A, Val(false)) # `Val(false)` turns off pivoting
    println("L+U:")
    display((L+U) .!= 0); println()
end



laplacian_1d(n) = (n+1)^2*Tridiagonal(
    fill( 1.0,n-1), # subdiagonal
    fill(-2.0,n),   # diagonal
    fill( 1.0,n-1)  # superdiagonal
)

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Δ,Id) + kron(Id,Δ)
end

"""
    nested_dissection(n) -> p

Compute the nested dissection permutation for `laplacian_2d(n)`.

Note that `n` must be of the form `n = 2^k - 1`.

For convenience, this function works with two-dimensional indices `i1,i2`
rather than the one-dimensional index `i = i1 + n*(i2-1)`. The output of this
function therefore needs to be converted using the
`matrix_to_vector_indices()` function provided below before it can be used to
permute a matrix.
"""
function nested_dissection(n)
    @assert isodd(n)
    n == 1 && return [(1,1)]

    nn = n÷2
    pp = nested_dissection(nn)
    return [
        [(i[1]     , i[2]     ) for i in pp];  # top left quadrant
        [(i[1]+nn+1, i[2]     ) for i in pp];  # bottom left quadrant
        [(i[1]     , i[2]+nn+1) for i in pp];  # top right quadrant
        [(i[1]+nn+1, i[2]+nn+1) for i in pp];  # bottom right quadrant
        [(nn+1,      i2) for i2 in 1:nn];      # left horizontal separator
        [(nn+1, nn+1+i2) for i2 in 1:nn];      # right horizontal separator
        [(i1, nn+1) for i1 in 1:n]             # vertical separator
    ]
end

matrix_to_vector_indices(n,p) = [i[1] + n*(i[2]-1) for i in p]

"""
    runtimes(n = 127)

Print the runtimes of LU factorisation applied to `P*laplacian_2d(n)*P` for
`P in [I, nested_dissection, AMD]`.

Note that `n` must be of the form `n = 2^k - 1`.
"""
function runtimes(n = 127)
    Δ = -laplacian_2d(n)

    t = @elapsed( cholesky(Δ; perm = 1:n^2) ) # cholesky == lu for symmetric matrices
    println("         Original: ", round(t, sigdigits=3), " sec")

    p = matrix_to_vector_indices(n,nested_dissection(n))
    t = @elapsed( cholesky(Δ; perm = p) )
    println("Nested dissection: ", round(t, sigdigits=3), " sec")

    t = @elapsed( cholesky(Δ) )
    println("              AMD: ", round(t, sigdigits=3), " sec")
end

"""
    sparsity_pattern(n, perm = "")

Plot the sparsity pattern of `L+U`, where `L`,`U` is the LU factorisation of
`Δ = laplacian_2d(n)`.

`perm in ["","nd","amd"]` denotes the permutation to apply to `Δ` before
the factorisation.

Note that for `perm == "nd"`, `n` must be of the form `n = 2^k - 1`.
"""
function sparsity_pattern(n,perm = "")
    A = -laplacian_2d(n)
    if perm == "nd"
        p = matrix_to_vector_indices(n, nested_dissection(n))
        A = A[p,p]
    end
    if perm == "amd"
        p = amd(A)
        A = A[p,p]
    end

    L,U = lu(Matrix(A), Val(false)) # `Val(false)` disables pivoting for stability

    clf()
    spy(L.+U, marker="s", c="r", ms=3e2*n^(-2))
    spy(  A , marker="s", c="k", ms=3e2*n^(-2))
    display(gcf())
end

