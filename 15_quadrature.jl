using PyPlot
using LaTeXStrings
using FastGaussQuadrature

function newton_cotes(n)
    x = LinRange(-1,1,n)
    w = quadweights(x)
    return x,w
end

function clenshaw_curtis(n)
    x = @. cos(π*(0:n-1)/(n-1))
    w = quadweights(x)
    return x,w
end

function convergence()
    x0 = 1.2
    f = x->1/(x+x0)
    I = log(1+x0) - log(-1+x0)      # exact integral of f on [-1,1]
    r = abs(x0 + sqrt(x0^2 - 1))  # rate of convergence
    n = 2:80

    clf()
    for (name,quadrule) in (
        ("NC", newton_cotes),
        ("CC", clenshaw_curtis),
        ("GL", gausslegendre),
    )
        errors = [begin
            x,w = quadrule(n)
            abs(sum(w.*f.(x)) - I)
        end for n in n]
        semilogy(n, errors, label=name)

        if quadrule == clenshaw_curtis
            nn = (10,45)
            semilogy(nn, r.^.-nn, "k--")
        end
        if quadrule == gausslegendre
            nn = (10,28)
            semilogy(nn, (r^2).^.-nn, "k-.")
        end
    end
    legend(loc="best")
end


"""
   legendre(x, n = length(x)) -> A

Assemble the matrix `A[i,j] = P_{j-1}(x[i])` with `P_j(x)` the `j`th Legendre polynomial.
"""
function legendre(x, n = length(x))
    n == 1 && return reshape( one.(x), (:,1) )
    n == 2 && return [ one.(x)  x ]
    A = [ one.(x)  x  zeros(length(x),n-2) ]
    for k = 2:n-1
        @. A[:,k+1] = ( (2k-1) * x * A[:,k] - (k-1) * A[:,k-1] ) / k
    end
    return A
end

"""
    quadweights(x) -> w

Determine quadrature weights `w` such that the quadrature rule `(x,w)` is exact
for all polynomials `p` with `degree(p) <= length(x)-1`.
"""
function quadweights(x)
    n = length(x)
    A = transpose(legendre(x))
    b = zeros(n); b[1] = 2
    return A \ b
end


function substitution_example()
    clf()

    subplot(1,2,1)
    title(L"\sqrt{1-x^2}")
    f = x->sqrt(1-x^2)
    n = round.(Int, 10.0.^LinRange(0,3,30))
    errors = [begin
        x,w = clenshaw_curtis(n)
        abs(sum(w.*f.(x)) - π/2)
    end for n in n]
    loglog(n, errors)
    loglog(n, 1e-1.*inv.(n).^3, "k--")

    subplot(1,2,2)
    title(L"\cos(\theta)^2")
    f = θ->cos(θ)^2
    n = 1:20
    errors = [begin
        x,w = clenshaw_curtis(n)
        θ = π/2*x
        wθ = π/2*w
        abs(sum(wθ.*f.(θ)) - π/2)
    end for n in n]
    semilogy(n, errors)

    tight_layout()
end



function midpoint(f,m)
    x = LinRange(-1,1,2m+1)[2:2:end]
    return 2/m*sum(f.(x))
end

function trapezoidal(f,m)
    x = LinRange(-1,1,m+1)
    return 1/m*(f(-1) + 2*sum(f.(x[2:end-1])) + f(1))
end

function simpson(f,m)
    x = LinRange(-1,1,2m+1)
    return 1/(3m)*(f(-1) + 4*sum(f.(x[2:2:end-1])) + 2*sum(f.(x[3:2:end-2])) + f(1))
end

function composite_convergence()
    x0 = 1.2
    f = x->1/(x+x0)
    I = log(1+x0) - log(-1+x0)      # exact integral of f on [-1,1]
    r = abs(x0 + sqrt(x0^2 - 1))  # rate of convergence
    m = round.(Int, 10.0.^LinRange(0,4, 30))

    clf()
    for (name,quadfun) in (
        ("Midpoint", midpoint),
        ("Trapezoidal", trapezoidal),
        ("Simpson", simpson),
    )
        errors = [begin
            abs(quadfun(f,m) - I)
        end for m in m]
        loglog(m, errors, label=name)

        mm = (1e1,1e4)
        loglog(mm, mm.^-2, "k--")
        loglog(mm, mm.^-4, "k-.")
    end
    legend(loc="best")
end
