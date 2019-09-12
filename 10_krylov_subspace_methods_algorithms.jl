using LinearAlgebra

function krylov_vectors(A,b,n)
    N = length(b)
    @assert size(A) == (N,N)

    V = zeros(N,n)
    for k = 1:n
        V[:,k] = b
        b = A*b
    end
    return V
end

function normalised_krylov_vectors(A,b,n)
    V = krylov_vectors(A,b,n)
    return V ./ [norm(V[:,k]) for k = 1:n]'
end

function gmres_unstable(A,b,n)
    V = krylov_vectors(A,b,n)
    y = ((A*V) \ b)
    return V * y
end

function arnoldi(A,b,n)
    N = length(b)
    @assert size(A) == (N,N)

    Q = zeros(N,n+1)
    H = zeros(n+1,n)

    Q[:,1] = b / norm(b)
    for k = 1:n
        q̃ = A*Q[:,k]
        for l = 1:k
            H[l,k] = dot(Q[:,l],q̃)
            q̃ = q̃ - H[l,k]*Q[:,l]
        end
        H[k+1,k] = norm(q̃)
        Q[:,k+1] = q̃ / H[k+1,k]
    end
    return Q,H
end

function gmres_slow(A,b,n)
    Q,H = arnoldi(A,b,n)
    y = (A*Q[:,1:end-1]) \ b
    return Q[:,1:end-1]*y
end

function gmres(A,b,n)
    Q,H = arnoldi(A,b,n)
    b̂ = zeros(n+1); b̂[1] = norm(b)
    y = H \ b̂
    #      ^ This does not exploit the Hessenberg structure of H.
    #        Use better algorithm in practice!
    return Q[:,1:end-1]*y
end

using Test, Random

function test(gmres_implementation)
    @testset for N = 2:20
        A = rand(N,N)
        x = rand(N)
        x̃ = gmres_implementation(A,A*x,N)
        @test x̃ ≈ x
    end
end

function test()
    Random.seed!(42)
    @testset "gmres_implementations" begin
        @testset "gmres_unstable" begin test(gmres_unstable); end
        @testset "gmres_slow" begin test(gmres_slow); end
        @testset "gmres" begin test(gmres); end
    end
end


function lanczos(A,b,n)
    N = length(b)
    @assert size(A) == (N,N)

    Q = zeros(N,n+1)
    H = zeros(n+1,n) # Use `Tridiagonal` in practice!

    Q[:,1] = b / norm(b)
    for k = 1:n
        q̃ = A*Q[:,k]
        H[k,k] = dot(Q[:,k],q̃)
        if k == 1
            q̃ = q̃ - H[k,k]*Q[:,k]
        else
            q̃ = q̃ - H[k,k]*Q[:,k] - H[k-1,k]*Q[:,k-1]
        end
        H[k+1,k] = norm(q̃)
        if k < n; H[k,k+1] = H[k+1,k]; end
        Q[:,k+1] = q̃ / H[k+1,k]
    end
    return Q,H
end

function minres(A,b,n)
    Q,H = lanczos(A,b,n)
    b̂ = zeros(n+1); b̂[1] = norm(b)
    y = H \ b̂   # Use better algorithm in practice
    return Q[:,1:end-1]*y
end

function test_minres()
    @testset "minres" begin
        @testset for N = 2:9
            A = rand(N,N)
            A += A' # Make sure A is symmetric
            x = rand(N)
            x̃ = minres(A,A*x,N)
            @test x̃ ≈ x
        end
    end
end
