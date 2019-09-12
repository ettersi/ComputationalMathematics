function laplacian_2d_tedious(n)
    i = Int[]
    j = Int[]
    v = Float64[]

    for k1 = 1:n, k2 = 1:n
        # Diagonal entry
        push!(i, k1 + n*(k2-1))
        push!(j, k1 + n*(k2-1))
        push!(v, -4.0)

        # Off-diagonal entries
        if k1 < n
            push!(i, k1+1 + n*(k2-1))
            push!(j, k1 + n*(k2-1))
            push!(v, 1)
        end

        if k1 > 1
            push!(i, k1-1 + n*(k2-1))
            push!(j, k1 + n*(k2-1))
            push!(v, 1)
        end

        if k2 < n
            push!(i, k1 + n*(k2+1-1))
            push!(j, k1 + n*(k2-1))
            push!(v, 1)
        end

        if k2 > 1
            push!(i, k1 + n*(k2-1-1))
            push!(j, k1 + n*(k2-1))
            push!(v, 1)
        end
    end

    return sparse(i,j,v)
end

laplacian_1d(n) = sparse(Tridiagonal(
    fill(-1.0,n-1), # subdiagonal
    fill( 2.0,n),   # diagonal
    fill(-1.0,n-1)  # superdiagonal
))

function laplacian_2d_clever(n)
    Δ = laplacian_1d(n)
    Id = sparse(I,n,n)
    return kron(Δ,Id) + kron(Id,Δ)
end
