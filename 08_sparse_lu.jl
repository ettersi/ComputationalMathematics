using PyPlot
using LinearAlgebra

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