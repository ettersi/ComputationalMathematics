# Run `] add FFTW` to install the FFTW package
using FFTW
using PyPlot
using LinearAlgebra
using SparseArrays

function solve_poisson(f)
    @assert size(f,1) == size(f,2)
    n = size(f,1)
    f̂ = FFTW.r2r(f, FFTW.RODFT00)
    λ = @.( (n+1)^2 * (2*cos(π*(1:n)/(n+1)) - 2) )
    û = f̂./(.- λ' .- λ)./(2*(n+1))^2
    u = FFTW.r2r(û, FFTW.RODFT00)
end

function example()
    n = 100
    x = LinRange(0,1, n+2)[2:end-1]
    f = (x,y) -> (0.3 < sqrt((x-0.5)^2 + (y-0.5)^2) < 0.4) * (abs(atan(x-0.5,y-0.5)) > π/4)
    u = solve_poisson(f.(x',x))

    clf()
    subplot(1,2,1)
    imshow(f.(x',x), origin="bottom left")
    subplot(1,2,2)
    imshow(u, origin="bottom left")
    colorbar()
end


##########
# Testing

using Test

function test()
    @testset for n = 2:5
        u = rand(n,n)
        Δ = laplacian_1d(n)
        f = -Δ*u - u*Δ
        ũ = solve_poisson(f)
        @test u ≈ ũ
    end
end

laplacian_1d(n) = sparse(Tridiagonal(
    fill((n+1)^2,n-1), # subdiagonal
    fill(-2*(n+1)^2,n),   # diagonal
    fill((n+1)^2,n-1)  # superdiagonal
))
