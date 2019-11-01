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
    x1,y1,dx1,dy1 = euler_step(x0,y0,dx0,dy0,dt/2)
    v1 = sqrt(dx1^2 + dy1^2)
    return (
        x0 + dx1*dt,
        y0 + dy1*dt,
        dx0 - v1*dx1*dt,
        dy0 - v1*dy1*dt - dt
    )
end

function heun_step(x0,y0,dx0,dy0,dt)
    x1,y1,dx1,dy1 = midpoint_step(x0,y0,dx0,dy0,dt*2/3)
    v0 = sqrt(dx0^2 + dy0^2)
    v1 = sqrt(dx1^2 + dy1^2)
    return (
        x0 + (dx0/4 + dx1*3/4)*dt,
        y0 + (dy0/4 + dy1*3/4)*dt,
        dx0 - (v0*dx0/4 + v1*dx1*3/4)*dt,
        dy0 - (v0*dy0/4 + v1*dy1*3/4)*dt - dt
    )
end

function semi_implicit_euler_step(x0,y0,dx0,dy0,dt)
    v0 = sqrt(dx0^2 + dy0^2)
    return x0+dx0*dt, y0+dy0*dt, dx0/(1+v0*dt), (dy0-dt)/(1+v0*dt)
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
        ("midpoint", midpoint_step),
        ("Heun", heun_step),
        ("semi-implicit Euler", semi_implicit_euler_step),
    )
        x,y,dx,dy = integrate(dx0,dy0,T,n,step)
        plot(x,y, label=name)
    end
    legend(loc="best")
    xlabel("x")
    ylabel("y")
end


function convergence()
    # Define the trajectory parameters
    T = 2
    dx0,dy0 = 1,1
    n = round.(Int,10.0.^LinRange(0,3,31))

    # Compute reference solution
    x,y,dx,dy = integrate(dx0,dy0,T,10000, heun_step)
    x̂,ŷ = x[end],y[end]
    dx̂,dŷ = dx[end],dy[end]

    # Compute and plot errors
    clf()
    for (name, step) in (
        ("Euler", euler_step),
        ("midpoint", midpoint_step),
        ("Heun", heun_step),
        ("semi-implicit Euler", semi_implicit_euler_step),
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
    xlabel("n")
    ylabel("error")
end

function stability()
    # Define the trajectory parameters
    dx0,dy0 = 0,-0.99
    n = 100

    # Compute the trajectories and plot |dy + 1|
    clf()
    for (name, step, dt) in (
        ("Euler", euler_step, 1),
        ("midpoint", midpoint_step, 1),
        ("Heun", heun_step, 2.5127453266183286240237/2),
        ("semi-implicit Euler", semi_implicit_euler_step, 1e3),
    )
        x,y,dx,dy = integrate(dx0,dy0,dt*n,n,step)
        semilogy(abs.(dy .+ 1), label=name)
    end
    legend(loc="best")
    xlabel("step number")
    ylabel("|dy + 1|")
end
