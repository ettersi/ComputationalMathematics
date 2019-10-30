using LinearAlgebra
using PyPlot

using PyCall
FuncAnimation = pyimport("matplotlib.animation").FuncAnimation

laplacian(n) = (n+1)^2*Tridiagonal(
    fill( 1,n-1), # subdiagonal
    fill(-2,n),   # diagonal
    fill( 1,n-1)  # superdiagonal
)

explicit_euler(A,u0,dt) = u0 .+ A*u0.*dt
implicit_euler(A,u0,dt) = (I - A*dt) \ u0
explicit_midpoint(A,u0,dt) = u0 .+ A*explicit_euler(A,u0,dt/2).*dt
implicit_midpoint(A,u0,dt) = u0 .+ A*implicit_euler(A,u0,dt/2).*dt

function example()
    # Initial conditions
    u0 = x -> abs(x-0.5) < 0.25

    # Simulation parameters
    n = 100
    step = explicit_euler
    dt = 0.1/(n+1)^2

    A = laplacian(n)
    x = LinRange(0,1,n+2)
    u = u0.(x[2:end-1])

    # Create animation
    simtime_to_realtime = 500
    steps_per_frame = round(Int, inv(10*dt*simtime_to_realtime), RoundUp)
    seconds_per_frame = simtime_to_realtime*dt*steps_per_frame

    clf()
    p, = plot(x, fill(NaN,n+2))
    xlim([-0.1,1.1])
    ylim([-0.1,1.1])

    t = 0.0
    FuncAnimation(
        gcf(),
        i->begin
            for i = 1:steps_per_frame
                u = step(A,u,dt)
            end
            p.set_ydata([0;u;0])
            return (p,)
        end,
        interval = seconds_per_frame*1000,
        blit=true,
        init_func=()->(p,)
    )
end

function convergence()
    u = (x,t)->sin(π*x)*exp(-π^2*t)
    m = [10,100,1000,10000,100000]
    n = round.(Int,sqrt.(m))
    step = implicit_euler
    T = 1

    errors = [begin
        x = LinRange(0,1,n+2)[2:end-1]
        A = laplacian(n)
        ũ = u.(x,0)
        for i = 1:m
            ũ = step(A,ũ,T/m)
        end
        norm(ũ .- u.(x,T),2)/sqrt(n+1)
    end for n in n, m in m]

    display(round.(log10.(errors), digits=1))
end
