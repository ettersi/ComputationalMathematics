"""
    mylu(A) -> L,U

Compute the unpivoted LU factorisation of `A`.
"""
function mylu(A)
    @assert size(A,1) == size(A,2) "Argument A in mylu(A) must be a square matrix"
    # TODO: your code here!
end

"""
    tril_solve(L,b) -> x

Compute `x = L\b` where `L` is lower-triangular.

Note that you may *not* assume `L[i,i] == 1` here.
"""
function tril_solve(L,b)
    @assert size(L,1) == size(L,2) "Argument L in tril_solve(L,b) must be a square matrix"
    @assert size(L,1) == length(b) "Sizes of arguments L,b in tril_solve(L,b) don't match"
    # TODO: your code here!
end

"""
    triu_solve(U,b) -> x

Compute `x = U\b` where `U` is upper-triangular.

The provided implementation transforms this linear system into an upper-triangular
one and then calls `tril_solve()`.
"""
triu_solve(U,b) = reverse(tril_solve(U[end:-1:1,end:-1:1], reverse(b)))

"""
    solve(A,b) -> x

Compute `x == A\b` using the `mylu()` factorisation.
"""
function solve(A,b)
    # TODO: your code here!
end


##########
# Testing

using Test, Random

function test()
    @testset "Sheet 2" begin
        test_LU()
        test_tril_solve()
        test_solve()
    end
end

function test_LU()
    @testset "mylu()" begin
        Random.seed!(42)
        @testset for n = 2:10
            for i = 1:5
                A = rand(n,n)
                L,U = mylu(A)

                @test all(diag(L) .≈ 1)
                @test all(triu(L,1) .== 0)
                @test all(tril(U,-1) .== 0)
                @test L*U ≈ A
            end
        end
    end
end

function test_tril_solve()
    @testset "tril_solve()" begin
        Random.seed!(42)
        @testset for n = 2:10
            for i = 1:5
                L = tril(rand(n,n))
                x = rand(n)

                @test x ≈ tril_solve(L,L*x)
            end
        end
    end
end

function test_solve()
    @testset "solve()" begin
        Random.seed!(42)
        @testset for n = 2:10
            for i = 1:5
                A = rand(n,n)
                x = rand(n)

                @test x ≈ solve(A,A*x)
            end
        end
    end
end


######################################
# Accuracy and performance assessment

using LinearAlgebra

function wilkinson_matrix(n)
    A = I + tril(fill(-1.0,(n,n)),-1)
    A[:,end] .= 1
    return A
end

# TODO: your code here!
#  - Compute the relative forward errors as described in Exercise 2.1.
#  - Measure the runtimes as described in Exercise 2.3.
