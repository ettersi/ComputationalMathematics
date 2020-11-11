using PyPlot
using Roots # Install by typing `] add Roots`

function simpson_step(f,y0,t)
    # TODO: Your code here!
end

function implicit_simpson_step(f,y0,t)
    # TODO: Your code here!
end

function integrate(f,y0,T,n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function convergence()
    f = y->y^2
    y0 = 1.0
    T = 0.5
    y = t-> y0/(1-y0*t)

    clf()
    n = round.(Int, 10.0.^LinRange(1,3,30))
    for (name,step) in (
        ("explicit", simpson_step),
        # ("implict", implicit_simpson_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n).^2, "k--")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end