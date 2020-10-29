using PyPlot

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function midpoint_step(f,y0,t)
    ỹ = y0 + f(y0)*t/2
    return y0 + f(ỹ)*t
end

function rk4_step(f,y0,t)
    f0 = t*f(y0)
    f1 = t*f(y0 + f0/2)
    f2 = t*f(y0 + f1/2)
    f3 = t*f(y0 + f2)
    return y0 + f0/6 + f1/3 + f2/3 + f3/6
end

function integrate(f,y0,T,n,step)
    y = Vector{float(typeof(y0))}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function example()
    f = y->[ y[2], -y[1] ]
    y0 = [ 1.0, 0.0 ]
    n = 1000
    t = LinRange(0,2π,n)

    clf()
    tt = LinRange(0,t[end],1000)
    plot(tt, cos.(tt), "k", label="exact")
    for (name,step) in (
        ("Euler", euler_step),
        # ("midpoint", midpoint_step),
        # ("RK4", rk4_step),
    )
        ỹ = integrate(f,y0,t[end],n, step)
        plot(t, [ỹ[i][1] for i = 1:n], label=name)
    end
    xlabel(L"t")
    ylabel(L"y(t)")
    legend(frameon=false)
    display(gcf())
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
        # ("midpoint", midpoint_step),
        # ("RK4", rk4_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n), "k--")
    # loglog(n, inv.(n).^2, "k-.")
    # loglog(n, inv.(n).^4, "k-.")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end

function nsteps()
    λ = 1.0
    f = y->λ*y
    y0 = one(λ)
    y = t->exp(λ*t)
    T = LinRange(0,3,11)
    τ = 1e-3

    n = [begin
        n = 2
        while true
            n = round(Int, n*1.3)
            ỹ = integrate(f,y0,T,n, euler_step)
            t = LinRange(0,T,n)
            if maximum(abs(y(t[i]) - ỹ[i]) for i = 1:n) < τ
                break
            end
        end
        n
    end for T in T]

    clf()
    plot(T,n)
    # semilogy(T,n)
    xlabel(L"T")
    ylabel(L"$n$ required to achieve error tolerance")
    display(gcf())
end
