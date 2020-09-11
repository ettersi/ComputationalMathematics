using PyPlot
using FFTW
using BernsteinEllipses
using FastGaussQuadrature

chebyshev_points(n) = cos.(LinRange(0,π,2n+1)[2:2:end-1])
function clenshaw_curtis(f,a,b,n)
    # WARNING:
    # This function uses some advanced techniques to evaluate the
    # Clenshaw-Curtis quadrature rule efficiently. You are not expected to
    # understand how this function works.

    ϕ = x -> a * (1-x)/2 + b * (1+x)/2
    x = chebyshev_points(n)
    c = FFTW.r2r(f.(ϕ.(x)),FFTW.REDFT10) ./ n
    return (
        c[1]
        +
        reduce(+, 2*c[k+1]/(1-k^2) for k = 2:n-1 if iseven(k); init=0.0)
    ) * (b-a)/2
end
function sample_chebyshev_interpolant(f,n,nn)
    # Input:
    #    f : function to interpolate
    #    n : number of interpolation points
    #   nn : number of evaluation points
    #
    # Output:
    #   xx : evaluation points of the interpolant
    #  pxx : values of the interpolant at xx

    # WARNING:
    # This function uses some advanced techniques to sample the
    # Chebyshev interpolants efficiently. You are not expected to
    # understand how this function works.

    x = chebyshev_points(n)
    c = FFTW.r2r(f.(x),FFTW.REDFT10)./(2n)

    xx = chebyshev_points(nn)
    pxx = FFTW.r2r([c; zeros(nn-n)], FFTW.REDFT01)

    return xx,pxx
end

function clenshaw_curtis_abs()
    n = round.(Int, 10.0.^LinRange(0,3,100))
    k = (1,3)
    c = ("C0", "C2")

    clf()
    for (k,c) in zip(k,c)
        a,b = -1,1
        f = x -> abs(x)^k
        I_ref = clenshaw_curtis(f, a,b, 2^14)
        errors = abs.(clenshaw_curtis.(f,a,b,n) .- I_ref)

        label = ifelse(k > 1, latexstring("|x|^$k"), L"|x|")
        loglog(n, errors, "$c", label=label)

        nn = (1e1,1e3)
        loglog(nn, 4e0.*float.(nn).^-(k+1), "$c--", label=latexstring("O(n^{-$(k+1)})"))
    end
    xlabel(L"n")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end

function chebyshev_error_abs()
    n = 100
    f = x->abs(x)
    xx,pxx = sample_chebyshev_interpolant(f,n,1024)

    clf()
    plot(xx, pxx .- f.(xx), "C0-")
    display(gcf())
end

function clenshaw_curtis_runge()
    n = 1:100
    a = (4,25)
    c = ("C0", "C2")
    s = (1e-2,1e-3)
    nn = ((20,60), (50,100))

    clf()
    for (a,c,s,nn) in zip(a,c,s,nn)
        f = x -> inv(1 + a*x^2)
        r = BernsteinEllipses.radius(im/sqrt(a))
        I_ref = clenshaw_curtis(f,-1,1,1024)
        errors = abs.(clenshaw_curtis.(f,-1,1,n) .- I_ref)
        semilogy(n, errors, "$c-", label=latexstring("\\frac{1}{1 + $a x^2}"))
        semilogy(nn, s.*r.^.-nn, "$c--", label=latexstring("O(|\\phi^{-1}($(1/sqrt(a)) i)|^{-n})"))
    end
    xlabel(L"n")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end

midpoint(f,a,b,m) = sum(f.(LinRange(a,b,2m+1)[2:2:end-1])) * (b-a)/ m
trapezoidal(f,a,b,m) = (f(a)/2 + sum(f.(LinRange(a,b,m+1)[2:end-1])) + f(b)/2) * (b-a)/m

function newton_cotes_convergence()
    a,b = 0,π
    f = x -> sin(x)
    Iref = clenshaw_curtis(f,a,b,1024)
    m = round.(Int, 10.0.^LinRange(1,3,50))

    clf()
    for (label,quad) = (
        ("midpoint", midpoint),
        ("trapezoidal", trapezoidal),
    )
        errors = abs.(quad.(f,a,b,m) .- Iref)
        loglog(m, errors, label=label)
    end
    mm = (4e1,1e3)
    loglog(mm, 4e0 .* mm.^(-2), "k--", label=L"O(m^{-2})")
    xlabel(L"m")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end

function gauss(f,a,b,n)
    # Quadrature rule for [-1,1]
    x,w = FastGaussQuadrature.gausslegendre(n)

    # Map to [a,b]
    x = (b+a)/2 .+ (b-a)/2 .* x
    w = (b-a)/2 .* w

    # Evaluate
    return sum(f.(x).*w)
end

function gauss_abs()
    n = round.(Int, 10.0.^LinRange(0,3,100))
    k = (1,3)
    c = ("C0", "C2")

    clf()
    for (k,c) in zip(k,c)
        a,b = -2,1
        f = x -> abs(x)^k
        I_ref = gauss(f, a,b, 2^14)
        errors = abs.(gauss.(f,a,b,n) .- I_ref)

        label = ifelse(k > 1, latexstring("|x|^$k"), L"|x|")
        loglog(n, errors, "$c", label=label)

        nn = (1e1,1e3)
        loglog(nn, 1e1.*float.(nn).^-(k+1), "$c--", label=latexstring("O(n^{-$(k+1)})"))
    end
    xlabel(L"n")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end

function gauss_runge()
    n = 1:100
    a = (4,25)
    c = ("C0", "C2")
    s = (2e1,1e1)
    nn = ((10,38), (10,90))

    clf()
    for (a,c,s,nn) in zip(a,c,s,nn)
        f = x -> inv(1 + a*x^2)
        r = BernsteinEllipses.radius(im/sqrt(a))
        I_ref = gauss(f,-1,1,1024)
        errors = abs.(gauss.(f,-1,1,n) .- I_ref)
        semilogy(n, errors, "$c-", label=latexstring("\\frac{1}{1 + $a x^2}"))
        semilogy(nn, @.(s*r^(-2*nn)), "$c--", label=latexstring("O(|\\phi^{-1}($(1/sqrt(a)) i)|^{-2n})"))
    end
    xlabel(L"n")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end

function clenshaw_curtis_vs_gauss()
    n = 1:100
    f = x -> inv(1+25x^2)
    I_ref = gauss(f,-1,1,1024)

    clf()
    for (label,quad) = (
        ("Clenshaw-Curtis", clenshaw_curtis),
        ("Gauss", gauss),
    )
        errors = abs.(quad.(f,-1,1,n) .- I_ref)
        semilogy(n, errors, label=label)
    end
    xlabel(L"n")
    ylabel("quadrature error")
    legend(frameon=false)
    display(gcf())
end
