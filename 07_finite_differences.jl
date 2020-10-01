using PyPlot
using LinearAlgebra

function laplacian_1d(n)
    return (n+1)^2 * Tridiagonal(
        fill( 1.0,n-1), # subdiagonal
        fill(-2.0,n),   # diagonal
        fill( 1.0,n-1)  # superdiagonal
    )
end

function solve_poisson_1d(f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian_1d(n)
    return x, -Δ\f.(x)
end

function example_1d()
    f = x -> π^2*sin(π*x)
    uref = x -> sin(π*x)
    x,u = solve_poisson_1d(f,4)

    clf()
    xx = LinRange(0,1,1000)
    plot(xx,uref.(xx), "k-", label="reference")
    plot([0;x;1],[0;u;0], label="FD")
    legend(frameon=false)
    xlabel(L"x")
    ylabel(L"u(x)")
    display(gcf())
end

function convergence_1d()
    # Define problem and solution
    f = x -> π^2 * sin(π*x)
    u = x -> sin(π*x)

    # Compute errors
    n = 2 .^ (1:15)
    error = [begin
        x,ũ = solve_poisson_1d(f,n)
        norm(ũ .- u.(x), 2)/sqrt(n+1)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end



using SparseArrays

function laplacian_2d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Id,Δ) + kron(Δ,Id)
end

function solve_poisson_2d(f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian_2d(n)
    b = vec(f.(x,x'))
    return x,reshape(-Δ\b, (n,n))
end

function example_2d()
    f = (x1,x2)->x1*x2
    x,u = solve_poisson_2d(f,300)

    clf()
    imshow(u, extent=(0,1,0,1), origin="bottom left")
    colorbar()
    display(gcf())
end

function convergence_2d()
    # Define problem and solution
    f = (x1,x2) -> 5*π^2 * sin(π*x1) * sin(2π*x2)
    u = (x1,x2) -> sin(π*x1) * sin(2π*x2)

    # Compute errors
    n = 2 .^ (1:9)
    error = [begin
        x,ũ = solve_poisson_2d(f,n)
        norm(ũ .- u.(x,x'), 2)/(n+1)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, 2e0*n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end



function laplacian_3d(n)
    Δ = sparse(laplacian_1d(n))
    Id = sparse(I,n,n)
    return kron(Id,Id,Δ) + kron(Id,Δ,Id) + kron(Δ,Id,Id)
end

function solve_poisson_3d(f, n)
    x = LinRange(0,1,n+2)[2:end-1]
    Δ = laplacian_3d(n)
    b = vec(f.(x,x',reshape(x,(1,1,n))))
    return x,reshape(-Δ\b, (n,n,n))
end

function convergence_3d()
    # Define problem and solution
    f = (x1,x2,x3) -> 3*π^2 * sin(π*x1) * sin(π*x2) * sin(π*x3)
    u = (x1,x2,x3) -> sin(π*x1) * sin(π*x2) * sin(π*x3)

    # Compute errors
    n = 2 .^ (1:5)
    error = [begin
        x,ũ = solve_poisson_3d(f,n)
        norm(ũ .- u.(x,x',reshape(x,(1,1,n))), 2)/(n+1)^(3/2)
    end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())
end