using PyPlot
using DifferentialEquations
using Random
using LinearAlgebra
using Statistics

using PyCall
FuncAnimation = pyimport("matplotlib.animation").FuncAnimation

function kinetic_energy(v)
    n = size(v,2)
    E = 0.0
    for i = 1:n
        E += (v[1,i]^2 + v[2,i]^2)/2
    end
    return E
end

function potential_energy(x,L)
    n = size(x,2)
    E = 0.0
    for i = 1:n
        for j = 1:i-1
            # Periodic boundary conditions:
            # When an leaves the simulation box on one side, we make it re-enter
            # on the other side. This means that in order to get consistent
            # forces, we must use the distance between the closest pair of
            # periodic images of atoms i and j.
            r2 = minimum(
                (x[1,i] - x[1,j] + sqrt(3)*L*s1)^2 +
                (x[2,i] - x[2,j] +         L*s2)^2
                for s1 = (-1,0,1), s2 = (-1,0,1)
            )
            ir2 = inv(r2)
            ir6 = ir2^3
            E += ir6^2 - 2*ir6
        end
    end
    return E
end

energy(x,v,L) = kinetic_energy(v) + potential_energy(x,L)


function initial_positions(n,L)
    x = Matrix{Float64}(undef, 2,2n^2)
    o = (L/2-n/2) .* (sqrt(3),1)
    for i = 1:n
        for j = 1:n
            x[:,2n*(i-1) + 2j-1] .= (sqrt(3)*(i-1  ), j-1  ) .+ o
            x[:,2n*(i-1) + 2j  ] .= (sqrt(3)*(i-0.5), j-0.5) .+ o
        end
    end
    return x
end

function initial_velocities(n,E)
    v = randn(2,2n^2)
    v .-= mean(v,dims=2)
    EE = kinetic_energy(v,)
    v .*= sqrt(E/EE)
    return v
end


function map_to_box!(x,L)
    # Impose periodic boundary conditions:
    # When an leaves the simulation box on one side, make it re-enter on the
    # other side.
    n = size(x,2)
    for i = 1:n
        if x[1,i] < 0
            x[1,i] += sqrt(3)*L
        elseif x[1,i] > sqrt(3)*L
            x[1,i] -= sqrt(3)*L
        end
        if x[2,i] < 0
            x[2,i] += L
        elseif x[2,i] > L
            x[2,i] -= L
        end
    end
    return x
end


function plot_positions(x,L)
    p, = plot(x[1,:],x[2,:], "o", ms = 10)
    xlim([0, sqrt(3)*L])
    ylim([0, L])
    gca().set_aspect("equal","box")
    return p
end


function assemble_problem(n,L,E,T)
    Random.seed!(42)
    x = initial_positions(n,L)
    v = initial_velocities(n,E)

    problem = HamiltonianProblem(
        (x,v,_)->energy(x,v,L),
        x,v, (0,T)
    )

    enforce_pbc = DiscreteCallback(
        (u,t,integrator) -> true,
        integrator -> begin
            x = integrator.u.x[1]
            map_to_box!(x,L)
        end
    )

    return problem, enforce_pbc
end


function example()
    # I could not figure out how to make the plot show up in VSCode, so you will
    # have to run this function from the Julia REPL.
    # You can do so using the following steps:
    #  - Run the Julia executable that comes with your Julia installation.
    #  - Type `cd("[path]")` where `[path]` to move to the directory of this file.
    #  - Type `include("11_molecular_dynamics.jl"); example()`
    # If all goes well, a new window should appear with the animation.
    # Feel free to contact me if you have issues getting this to work.

    n = 3
    L = 5*n
    E = 2n^2 * 2.0   # Energy of the system
                     # Increase this parameter to go from solid to liquid to gas
    dt = 0.01
    dt_per_frame = 0.1
    frames_per_second = 20

    problem, enforce_pbc = assemble_problem(n,L,E,Inf)

    integrator = init(
        problem,
        Midpoint(),
        adaptive=false, dt = dt,
        alias_u0 = true,
        save_on = false,
        callback=enforce_pbc
    )

    clf()
    p = plot_positions(integrator.u.x[1],L)
    FuncAnimation(
        gcf(),
        i->begin
            step!(integrator, dt_per_frame)
            p.set_data(integrator.u.x[1])
            return (p,)
        end,
        interval = 1000/frames_per_second,
        blit=true,
        init_func=()->(p,)
    )
end

function energy_conservation()
    n = 3
    L = 5*n
    E = 2n^2 * 0.1
    T = 10.0

    dt = 0.01
    method = Midpoint()         # Energy drifts
    # method = VerletLeapfrog()   # Energy oscillates but does not drift
    # dt = 0.06                   # Still stable for VerletLeapfrog()
    # dt = 0.07                   # Exponential blow-up for VerletLeapfrog()

    problem,enforce_pbc = assemble_problem(n,L,E,T)

    energy_values = SavedValues(Float64, NTuple{2,Float64})
    save_energies = SavingCallback(
        (u,t,integrator) -> (
            kinetic_energy(u.x[2]),
            potential_energy(u.x[1],L)
        ),
        energy_values
    )

    solve(
        problem,
        method,
        adaptive=false, dt = dt,
        alias_u0 = true,
        save_on = false,
        callback = CallbackSet(enforce_pbc,save_energies)
    )

    t = energy_values.t
    Ekin = [energy_values.saveval[i][1] for i = 1:length(t)]
    Epot = [energy_values.saveval[i][2] for i = 1:length(t)]

    clf()
    plot(t, Ekin.+Epot, label="total")

    # ---------------------------------------
    # Comment these lines to see the variation in the total energy more clearly
    plot(t, Ekin, label="kinetic")
    plot(t, Epot, label="potential")
    # ---------------------------------------

    legend(loc="best")
    xlabel("time")
    ylabel("energy")
    display(gcf())
end
