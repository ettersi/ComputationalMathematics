function classical_gram_schmidt_qr(A)
    m,n = size(A)
    Q = zeros(m,n)
    R = zeros(n,n)

    for k = 1:n
        R[1:k-1,k] = Q[:,1:k-1]' * A[:,k]
        q̃ = A[:,k] - Q[:,1:k-1] * R[1:k-1,k]
        R[k,k] = norm(q̃)
        Q[:,k] = q̃/R[k,k]
    end

    return Q,R
end

function modified_gram_schmidt_qr(A)
    m,n = size(A)
    Q = zeros(m,n)
    R = zeros(n,n)

    for k = 1:n
        q̃ = A[:,k]
        for l = 1:k-1
            R[l,k] = dot(Q[:,l], q̃)
            q̃ = q̃ - Q[:,l] * R[l,k]
        end
        R[k,k] = norm(q̃)
        Q[:,k] = q̃/R[k,k]
    end

    return Q,R
end

function householder_qr(A)
    m,n = size(A)
    Q = Matrix{Float64}(I, (m,m))
    R = copy(A)

    for k = 1:min(m,n)
        u = R[k:m,k]
        u[1] = u[1] + sign(u[1])*norm(u)
        u /= norm(u)
        Q[:,k:m] = Q[:,k:m] - 2*Q[:,k:m]*u*(u')
        R[k:m,:] = R[k:m,:] - 2*u*(u'*R[k:m,:])
    end

    return Q,R
end


using Test
using Random

function test_qr(qr,A)
    Q,R = qr(A)
    @test Q'*Q ≈ I
    @test norm(tril(R,-1)) <= sqrt(eps()) * norm(triu(R))
    @test Q*R ≈ A
end

function test_qr(qr)
    Random.seed!(42)
    test_qr(qr, rand(4,1))
    test_qr(qr, rand(4,2))
    test_qr(qr, rand(4,3))
    test_qr(qr, rand(4,4))
    test_qr(qr, rand(10,4))
end

function test()
    @testset "Classical Gram-Schmidt QR" begin test_qr(classical_gram_schmidt_qr); end
    @testset "Modified Gram-Schmidt QR" begin test_qr(modified_gram_schmidt_qr); end
    @testset "Householder QR" begin test_qr(householder_qr); end
end
