using PyPlot

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function midpoint_step(f,y0,t)
    ỹ = y0 + f(y0)*t/2
    return y0 + f(ỹ)*t
end

function trapezoidal_step(f,y0,t)
    ỹ = y0 + f(y0)*t
    return y0 + (f(y0) + f(ỹ))*t/2
end

function integrate(f,y0,T,n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function example()
    f = y->[ y[2], -y[1] ]
    y0 = [ 1.0, 0.0 ]
    n = 20
    t = LinRange(0,2π,n)

    clf()
    tt = LinRange(0,t[end],10000)
    plot(tt, cos.(tt), "k", label="ref")
    for (name,step) in (
        ("Euler", euler_step),
        ("midpoint", midpoint_step),
        ("trapezoidal", trapezoidal_step),
    )
        ỹ = integrate(f,y0,t[end],n, step)
        plot(t, [ỹ[i][1] for i = 1:n], label=name)
    end
    legend(loc="best")
end

function convergence()
    f = y->y^2
    y0 = 1.0
    T = 0.5
    y = t-> y0/(1-y0*t)

    clf()
    n = round.(Int, 10.0.^LinRange(0,3,30))
    for (name,step) in (
        ("Euler", euler_step),
        ("midpoint", midpoint_step),
        ("trapezoidal", trapezoidal_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n), "k--")
    loglog(n, inv.(n).^2, "k-.")
    legend(loc="best")
end
