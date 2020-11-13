using PyPlot
using LinearAlgebra
using Roots

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function trapezoidal_step(f,y0,t)
    f1 = t*f(y0)
    f2 = t*f(y0 + f1)
    return y0 + (f1 + f2)/2
end

function rk4_step(f,y0,t)
    f1 = t*f(y0)
    f2 = t*f(y0 + f1/2)
    f3 = t*f(y0 + f2/2)
    f4 = t*f(y0 + f3)
    return y0 + f1/6 + f2/3 + f3/3 + f4/6
end

function propagate(f,y0,T,n,step)
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
    tt = LinRange(0,t[end],1000)
    plot(tt, cos.(tt), "k", label="exact")
    for (name,step) in (
        ("Euler", euler_step),
        # ("trapezoidal", trapezoidal_step),
        # ("RK4", rk4_step),
    )
        ỹ = propagate(f,y0,t[end],n, step)
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
    n = round.(Int, 10.0.^LinRange(1,3,30))
    for (i,(name,step,p)) in enumerate((
        ("Euler", euler_step,1),
        # ("trapezoidal", trapezoidal_step,2),
        # ("RK4", rk4_step,4),
        # ("implicit Euler", implicit_euler_step,1),
        # ("implicit trapezoidal", implicit_trapezoidal_step,2),
    ))
        error = [begin
            ỹ = propagate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
        nn = (1e1,1e3)
        loglog(nn, inv.(nn).^p, "C$(i-1)--", label=latexstring("O(n^{-$p})"))
    end
    legend(frameon=false)
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
            ỹ = propagate(f,y0,T,n, euler_step)
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



function embedded_ET_step(f,y0,t)
    f1 = t*f(y0)
    y_euler = y0 + f1
    f2 = t*f(y_euler)
    y_trapezoidal = y0 + (f1+f2)/2
    return y_euler, y_trapezoidal
    # Use `norm(y_euler - y_trapezoidal)` as an approximation to the local error.
end

function propagate_adaptively(f,y0,T,τ,step,p)
    # p = order of consistency of step()[1]
    t = Vector{Float64}()
    y = Vector{typeof(y0)}()
    push!(t, 0.0)
    push!(y, y0)

    n_rejected = 0
    Δt = T
    while t[end] < T
        y1,y2 = step(f,y[end],Δt)
        q = (τ / norm(y1 - y2))^(1/p)
        if q >= 1
            push!(t,t[end]+Δt)
            push!(y,y2)
        else
            n_rejected += 1
        end
        Δt = min(0.9*q*Δt, T-t[end])
    end

    return t,y, n_rejected
end

function adaptive_rk_example_1()
    f = y->cos(y)^2
    y0 = -1.56
    T = 200
    τ = 1e-3

    t,y,n_rejected = propagate_adaptively(f,y0,T,τ, embedded_ET_step, 2)

    Δt = minimum(diff(t))
    n_fixed = round(Int, T/Δt)
    println("  # adaptive steps: ", length(t)-1)
    println("# fixed-size steps: ", n_fixed)
    println("             Ratio: ", round(n_fixed/(length(t)-1), sigdigits=3))
    println("  # rejected steps: ", n_rejected)

    clf()
    plot(t, zero.(t), "ko", ms=3, label=L"t_k")
    plot(t,y, "-", label=L"y(t)")
    xlabel(L"t")
    legend(frameon=false)
    display(gcf())
end

function adaptive_rk_example_2()
    f = y->-y
    y0 = 1.0
    T = 50
    τ = 1e-4

    if (explicit = true)
        t,y,_ = propagate_adaptively(f,y0,T,τ, embedded_ET_step, 2)
    else
        t,y,_ = propagate_adaptively(f,y0,T,τ, embedded_implicit_ET_step, 2)
    end

    clf()
    plot(t, fill(0.5, length(t)), "ko", ms=3, label=L"t_k")
    if (logy = false)
        semilogy(t,abs.(y), label=L"|y(t)|")
    else
        plot(t,y, label=L"y(t)")
    end
    xlabel(L"t")
    legend(frameon=false)
    display(gcf())
end

function stepsize()
    f = y->-y
    y0 = 1.0
    T = 50
    τ = 1e-6

    if (explicit = true)
        t,y,_ = propagate_adaptively(f,y0,T,τ, embedded_ET_step, 2)
    else
        t,y,_ = propagate_adaptively(f,y0,T,τ, embedded_implicit_ET_step, 2)
    end

    clf()
    if (show_ref = false)
        semilogy([t[1],t[end]],[2,2], "k", lw=0.5)
    end
    semilogy(t[2:end],diff(t))
    xlabel(L"t")
    ylabel(L"t_k - t_{k-1}")
    display(gcf())
end

function implicit_euler_step(f,y0,t)
    return find_zero(
        y -> y0 + f(y)*t - y,
        euler_step(f,y0,t)
    )
end

function implicit_trapezoidal_step(f,y0,t)
    f1 = t*f(y0)
    return find_zero(
        y -> y0 + (f1 + f(y)*t)/2 - y,
        trapezoidal_step(f,y0,t)
    )
end

function embedded_implicit_ET_step(f,y0,t)
    f1 = t*f(y0)
    y_trapezoidal = find_zero(
        y -> y0 + (f1 + f(y)*t)/2 - y,
        trapezoidal_step(f,y0,t)
    )
    y_quasi_euler = y0 + f(y_trapezoidal)*t
    return y_quasi_euler, y_trapezoidal
end


function stability_example()
    f = y->-y
    y0 = 1.0
    if (explicit = true)
        T = 10
        Δt = 1.2:0.2:2.2
        step = euler_step
        # step = trapezoidal_step
    else
        T = 50
        Δt = 1:5
        step = implicit_euler_step
        # step = implicit_trapezoidal_step
    end

    clf()
    t = LinRange(0,T,1000)
    plot(t, exp.(-t), "k", label=L"y(t)")
    for Δt = Δt
        n = round(Int,T/Δt)
        t = Δt.*(0:n)
        y = propagate(f,y0,n*Δt,n+1,step)
        plot(t,y, label=latexstring("\\Delta t = $Δt"))
    end
    xlabel(L"t")
    legend(frameon=false, loc="center left", bbox_to_anchor=(1,0.5))
    display(gcf())
end
