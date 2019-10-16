"""
    adaptive_simpson(f,(a,b),ε)

Compute the integral of `f(x)` over `[a,b]` up to an error `ε` using Simpson's
rule with the midpoint rule for error estimation.
"""
function adaptive_simpson(f,(a,b),ε)
    m = (a+b)/2
    return adaptive_simpson_recursion(f, (f(a),f(m),f(b)), (a,m,b), ε)
end

function adaptive_simpson_recursion(f,(fa,fm,fb),(a,m,b),ε)
    midpoint = (b-a)*fm
    simpson = (b-a)/6 * (fa + 4*fm + fb)

    if abs(midpoint-simpson) > ε
        ml = (a+m)/2
        mr = (m+b)/2
        return (
            adaptive_simpson_recursion(f, (fa, f(ml), fm), (a,ml,m), ε/2)
            +
            adaptive_simpson_recursion(f, (fm, f(mr), fb), (m,mr,b), ε/2)
        )
    else
        return simpson
    end
end

function example()
    # Function to integrate:
    f = x->1/(1 + 100*(x-0.1)^2)

    # Vector to store the evaluation points
    qx = Vector{Float64}()
    ff = x->( push!(qx,x); f(x) )

    # Perform the adaptive quadrature
    Q = adaptive_simpson(ff, (-1,1), 1e-2)

    # Plot function and evaluation points
    x = LinRange(-1,1,1000)
    clf()
    plot(x, f.(x), "k")
    plot(qx, f.(qx), "ko", ms=5)
end

function convergence()
    # Integrand and its anti-derivative
    f = x->1/(1 + 100*(x-0.1)^2)
    F = x->atan(10*(x-0.1))/10
    I = F(1) - F(-1)

    # Compute errors
    ε = 10.0.^LinRange(-8,0,100)
    errors = zeros(length(ε))
    nevals = zeros(Int, length(ε))
    for i = 1:length(ε)
        ff = x->( nevals[i] += 1; f(x) )
        Q = adaptive_simpson(ff,(-1,1),ε[i])
        errors[i] = abs(I-Q)
    end

    # Plot results
    clf()
    subplot(1,2,1)
    title("Error")
    loglog(ε, errors)
    loglog(ε,ε,"k")
    xlabel(L"\varepsilon")
    ylabel("error")

    subplot(1,2,2)
    title("Number of function evaluations")
    loglog(ε, nevals)
    ee = (ε[1],1e-2)
    loglog(ee,ee.^(-1/2), "k--", label=L"O(\varepsilon^{-1/2})")
    xlabel(L"\varepsilon")
    ylabel("# evaluations")
    legend(loc="best", frameon=false)

    gcf().set_size_inches(8,4)
    tight_layout()
end
