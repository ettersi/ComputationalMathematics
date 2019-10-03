using PyPlot

"""
    lagrange(xx,f, x)

Evaluate `p(x)` with `p` the polynomial interpolant through the data points `(xx[i],f[i])`.
"""
function lagrange(xx,f, x)
    @assert length(xx) == length(f)

    n = length(f)
    return sum(
        f[j] *
        prod(x .- xx[setdiff(1:n,j)]) /
        prod(xx[j] .- xx[setdiff(1:n,j)])
        for j = 1:n
    )
end

function example()
    xx = LinRange(-1,1,5)
    f = [1,0,2,1,2]
    p = x->lagrange(xx,f,x)

    x = LinRange(-1,1,1000)
    clf()
    plot(x, p.(x))
    plot(xx, f, "ko", ms=6)
end
