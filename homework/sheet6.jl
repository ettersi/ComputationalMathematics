using PyPlot

function euler_step(x0,y0,dx0,dy0,dt)
    v0 = sqrt(dx0^2 + dy0^2)
    return (
        x0 + dx0*dt,
        y0 + dy0*dt,
        dx0 - v0*dx0*dt,
        dy0 - v0*dy0*dt - dt
    )
end

function midpoint_step(x0,y0,dx0,dy0,dt)
    # TODO: your code here!
end

function heun_step(x0,y0,dx0,dy0,dt)
    # TODO: your code here!
end

function semi_implicit_euler_step(x0,y0,dx0,dy0,dt)
    # TODO: your code here!
end

function integrate(dx0,dy0,T,n,step)
    x = zeros(n+1)
    y = zeros(n+1)
    dx = zeros(n+1); dx[1] = dx0
    dy = zeros(n+1); dy[1] = dy0
    for i = 1:n
        x[i+1],y[i+1],dx[i+1],dy[i+1] = step(x[i],y[i],dx[i],dy[i],T/n)
    end
    return x,y,dx,dy
end

function trajectory()
    # Define the trajectory parameters
    T = 2
    dx0,dy0 = 1,1
    n = 20

    # Plot the solution trajectories
    clf()
    for (name, step) in (
        ("Euler", euler_step),
        # TODO: uncomment these lines once you have implemented the functions
        # ("midpoint", midpoint_step),
        # ("Heun", heun_step),
        # ("semi-implicit Euler", semi_implicit_euler_step),
    )
        x,y,dx,dy = integrate(dx0,dy0,T,n,step)
        plot(x,y, label=name)
    end
    legend(loc="best")
end


function convergence()
    # Define the trajectory parameters
    T = 2
    dx0,dy0 = 1,1
    n = round.(Int,10.0.^LinRange(0,3,31))

    # Compute reference solution
    # TODO: use a higher-order method to compute the reference solution
    x,y,dx,dy = integrate(dx0,dy0,T,10000, euler_step)
    x̂,ŷ = x[end],y[end]
    dx̂,dŷ = dx[end],dy[end]

    # Compute and plot errors
    clf()
    for (name, step) in (
        ("Euler", euler_step),
        # TODO: uncomment these lines once you have implemented the functions
        # ("midpoint", midpoint_step),
        # ("Heun", heun_step),
        # ("semi-implicit Euler", semi_implicit_euler_step),
    )
        error = [begin
            x,y,dx,dy = integrate(dx0,dy0,T,n,step)
            sqrt((x̂-x[end])^2 + (ŷ-y[end])^2)
        end for n in n]
        loglog(n,error, label=name)
    end

    # Plot reference lines
    nn = (1e1,1e3)
    loglog(nn, 2.0.*inv.(nn), "k--")
    loglog(nn, 0.5.*inv.(nn).^2, "k--")
    loglog(nn, 0.5.*inv.(nn).^3, "k--")

    legend(loc="best")
end

function stability()
    # Define the trajectory parameters
    dx0,dy0 = 0,-0.99
    n = 100

    # Compute the trajectories and plot |dy + 1|
    clf()
    for (name, step, dt) in (
        # TODO: uncomment these lines once you have implemented the functions
        # ("Euler", euler_step, TODO),
        # ("midpoint", midpoint_step, TODO),
        # ("Heun", heun_step, TODO),
        # ("semi-implicit Euler", semi_implicit_euler_step, 1e3),
    )
        x,y,dx,dy = integrate(dx0,dy0,dt*n,n,step)
        semilogy(abs.(dy .+ 1), label=name)
    end
    legend(loc="best")
end
