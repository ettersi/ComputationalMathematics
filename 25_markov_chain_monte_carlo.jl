using PyPlot
using FastGaussQuadrature

gfun(x1,x2) = exp(-10*(x1^2+x2^2-1)^2)
hfun(x1,x2) = (x1^2+x2)

function generate_samples(n,Δx)
    x1 = zeros(n)
    x2 = zeros(n)
    for i = 1:n-1
        x̃1 = x1[i] + Δx*randn()
        x̃2 = x2[i] + Δx*randn()
        R = gfun(x̃1,x̃2) / gfun(x1[i],x2[i])
        if rand() <= R
            x1[i+1] = x̃1
            x2[i+1] = x̃2
        else
            x1[i+1] = x1[i]
            x2[i+1] = x2[i]
        end
    end
    return x1,x2
end

function plot_samples()
    x = LinRange(-2,2,1000)
    clf()
    imshow(
        gfun.(x',x);
        extent=(x[1],x[end],x[1],x[end])
    )
    x1,x2 = generate_samples(1000,0.2)
    plot(x1,x2, "k")
end

function reference_value()
    x,w = gausslegendre(1000)
    x .*= 2; w .*= 2
    return sum(hfun(x1,x2)*gfun(x1,x2)*w1*w2 for (x1,w1) in zip(x,w), (x2,w2) in zip(x,w)) /
           sum(gfun(x1,x2)*w1*w2 for (x1,w1) in zip(x,w), (x2,w2) in zip(x,w))
end

function monte_carlo_value(Δx)
    N = 1_000_000
    x1,x2 = generate_samples(N,Δx)
    return sum(hfun.(x1,x2))/N
end

function step_length()
    Δx = 10.0.^LinRange(-3,2,101)
    I = reference_value()
    Q = monte_carlo_value.(Δx)

    clf()
    loglog(Δx, abs.(Q.-I)./I);
    xlabel("Δx")
    ylabel("relative error")
end
