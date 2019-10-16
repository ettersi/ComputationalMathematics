using PyPlot
using Roots

function explicit_euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function explicit_midpoint_step(f,y0,t)
    ỹ = y0 + f(y0)*t/2
    return y0 + f(ỹ)*t
end

function implicit_euler_step(f,y0,t)
    return find_zero(y->y0 + f(y)*t - y, (0,1))
    # `find_zero(f,(a,b))` determines the root of `f` in the interval `(a,b)`.
    # The choice `(a,b) = (0,1)` is problem-specific and will not work in general.
end

function implicit_midpoint_step(f,y0,t)
    ỹ = find_zero(y->y0 + f(y)*t/2 - y, (-1,1))
    # `find_zero(f,(a,b))` determines the root of `f` in the interval `(a,b)`.
    # The choice `(a,b) = (-1,1)` is problem-specific and will not work in general.
    return y0 + f(ỹ)*t
end

function integrate(f,y0,T,m,step)
    y = Vector{typeof(y0)}(undef,m)
    y[1] = y0
    for i = 2:m
        y[i] = step(f,y[i-1],T/(m-1))
    end
    return y
end

function example_explicit()
    f = y -> -y
    y0 = 1.0
    T = 10

    clf()
    t = LinRange(0,T,1000)
    plot(t, exp.(.-t), "k", label="ref")
    for m in (4,5,6,7,8)
        t = LinRange(0,T,m)
        ỹ = integrate(f,y0,T,m, explicit_midpoint_step)
        plot(t, [ỹ[i][1] for i = 1:m], label="m = $m")
    end
    ylims = ylim()
    @show ylims
    ylim(clamp.(ylims, -2.4,2.4))
    legend(loc="best")
end

function harmonic_oscillator()
    f = y -> [y[2],-y[1]]
    y0 = [1.0,0.0]
    m = 30
    T = 10π
    t = LinRange(0,T,m)

    clf()
    for (name,step) in (
        ("Euler", explicit_euler_step),
        ("midpoint", explicit_midpoint_step),
    )
        ỹ = integrate(f,y0,t[end],m, step)
        semilogy(t, [abs(ỹ[i][1]) for i = 1:m], label=name)
    end
    z = im*T/(m-1)
    semilogy(t, abs(1+z).^(0:m-1), "k--")
    semilogy(t, abs(1+z+z^2/2).^(0:m-1), "k--")
    legend(loc="best")
end

function example_implicit()
    f = y -> -y
    y0 = 1.0
    T = 10

    clf()
    t = LinRange(0,T,1000)
    plot(t, exp.(.-t), "k", label="ref")
    for m in (4,6,8)
        t = LinRange(0,T,m)
        ỹ = integrate(f,y0,T,m, implicit_midpoint_step)
        plot(t, [ỹ[i][1] for i = 1:m], label="m = $m")
    end
    legend(loc="best")
end

function convergence()
    f = y->-y
    y0 = 1.0
    T = 10
    y = t->exp(-t)

    clf()
    m = round.(Int, 10.0.^LinRange(0,4,100))
    for (name,step) in (
        ("explicit Euler", explicit_euler_step),
        ("explicit midpoint", explicit_midpoint_step),
        ("implicit Euler", implicit_euler_step),
        ("implicit midpoint", implicit_midpoint_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,m, step)
            abs(y(T) - ỹ[end])
        end for m in m]
        loglog(m, error, label=name)
    end
    nn = (1e2,1e4)
    loglog(nn, 1e-2.*inv.(nn), "k--")
    loglog(nn, 4e-2.*inv.(nn).^2, "k-.")
    legend(loc="best")
end
