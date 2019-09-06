"""
    mylu(A) -> L,U

Compute the unpivoted LU factorisation of `A`.
"""
function mylu(A)
    @assert size(A,1) == size(A,2) "Argument A in mylu(A) must be a square matrix"
    n = size(A,1)
    L = Matrix{Float64}(I,(n,n))
    U = copy(A)
    for k = 1:n
        for i = k+1:n
            L[i,k] = U[i,k]/U[k,k]
            U[i,k] = 0
            for j = k+1:n
                U[i,j] = U[i,j] - L[i,k] * U[k,j]
            end
        end
    end
    return L,U
end

"""
    tril_solve(L,b) -> x

Compute `x = L\b` where `L` is lower-triangular.

Note that you may *not* assume `L[i,i] == 1` here.
"""
function tril_solve(L,b)
    @assert size(L,1) == size(L,2) "Argument L in tril_solve(L,b) must be a square matrix"
    @assert size(L,1) == length(b) "Sizes of arguments L,b in tril_solve(L,b) don't match"
    x = zeros(length(b))
    for i = 1:length(b)
        x[i] = (b[i] - dot(L[i,1:i-1],x[1:i-1]))/L[i,i]
    end
    return x
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
    L,U = mylu(A)
    y = tril_solve(L,b)
    return triu_solve(U,y)
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

relative_error(x,xx) = norm(x-xx,Inf) / norm(x)

function print_accuracies(A,x)
    println("mylu: ", relative_error(x,solve(A,A*x)))
    println("  lu: ", relative_error(x,lu(A)\(A*x)))
    println("  qr: ", relative_error(x,qr(A)\(A*x)))
end

function accuracy(n)
    println("Random matrix:")
    Random.seed!(42)
    A = rand(n,n)
    x = rand(n)
    print_accuracies(A,x)
    println()

    println("Wilkinson matrix:")
    Random.seed!(42)
    A = wilkinson_mat(n)
    x = rand(n)
    print_accuracies(A,x)
    println()
end

function timings(n)
    A = rand(n,n)
    println("mylu: ", @elapsed(mylu(A))," sec")
    println("  lu: ", @elapsed(  lu(A))," sec")
    println("  qr: ", @elapsed(  qr(A))," sec")
end
