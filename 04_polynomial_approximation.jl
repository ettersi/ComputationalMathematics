using PyPlot
using LinearAlgebra
using FFTW
using BernsteinEllipses  # Written by yours truly!

function lagrange(x,j,xx)
    #  x : interpolation points
    # xx : evaluation point
    return prod(
        ( xx - x[i] ) / ( x[j] - x[i] )
        for i = 1:length(x) if i != j
    )
end
function interpolate(x,f,xx)
    #  x : interpolation points
    #  f : interpolation values
    # xx : evaluation point
    return sum(f[k] * lagrange(x,k,xx) for k = 1:length(x))
end

function lagrange_example()
    x = [1,2,3,4]                # interpolation points
    f = rand(1:4,4)              # interpolation values
    xx = LinRange(0.5,4.5,200)   # evaluation points for plotting

    clf()
    plot(xx, interpolate.(Ref(x),Ref(f),xx))
    plot(x,f, "ko", ms = 4)
    display(gcf())
end

function piecewise_interpolation_example()
    a,b = 0,2π
    f = sin
    # f = x -> abs(x - sqrt(2))
    m = 2  # number of intervals
    n = 2  # polynomial degree
    y = LinRange(a,b,m+1)

    clf()
    for k = 1:m
        x = LinRange(y[k],y[k+1],n+1)
        p = xx -> interpolate(x,f.(x),xx)
        xx = LinRange(y[k],y[k+1],100)
        plot(xx, f.(xx), "k-", lw=0.5)
        plot(xx, p.(xx), lw=2)
        plot(x, f.(x), "ko", ms=4)
    end
    display(gcf())
end

function adaptive_piecewise_interpolation()
    a,b = 0,2π
    f = x -> abs(x - sqrt(2))
    n = 1   # polynomial degree
    τ = 0.2 # error tolerance

    clf()
    y, dy = a,b-a
    while y < b
        while true
            x = LinRange(y,y+dy,n+1)
            p = xx -> interpolate(x,f.(x),xx)
            xx = LinRange(y,y+dy,100)

            if norm(f.(xx) .- p.(xx),Inf) < τ
                plot(xx, f.(xx), "k-", lw=0.5)
                plot(xx, p.(xx), lw=2)
                plot(x, f.(x), "ko", ms=4)

                y += dy
                dy = min(2*dy, b-y)
                break
            end
            dy /= 2
        end
    end
    display(gcf())
end

function piecewise_interpolation_convergence()
    if true
        a,b = 0,2π
        f = sin

        c = [1e1, 5e0, 2e0, 1e0]
        plot_bigO = n->true
    else
        a,b = (-1,1)
        f = x->sign(x)*x^3

        c = [5e0, 8e-1]
        plot_bigO = n->n<=2
    end

    n = 1:4
    m = round.(Int, 10.0.^LinRange(1,3,20))

    clf()
    for (j,n) = enumerate(n)
        errors = zeros(length(m))
        for (i,m) = enumerate(m)
            y = LinRange(a,b,m+1)
            for k = 1:m
                x = LinRange(y[k],y[k+1],n+1)
                p = xx -> interpolate(x,f.(x),xx)

                xx = LinRange(y[k],y[k+1],100)
                errors[i] = max(errors[i], norm(f.(xx) .- p.(xx), Inf))
            end
        end
        loglog(m, errors, "C$(j-1)", label="n = $n")

        mm = [3e1,1e3]
        if plot_bigO(n)
            loglog(mm, c[j].*float.(mm).^(-n-1), "C$(j-1)--")
        end
    end
    xlabel(L"m")
    ylabel(L"\|f - p\|_{[a,b]}")
    legend(loc="upper right")
    display(gcf())
end

function node_polynomial()
    n = 5
    # n = 15

    clf()
    for (label,x) = (
        ( "uniform", LinRange(-1,1,2n+1)[2:2:end-1] ),
        # ( "Chebyshev", cos.(LinRange(0,π,2n+1)[2:2:end-1]) ),
    )
        l = xx->prod(xx .- x)
        xx = LinRange(-1,1,1000)
        if true
            plot(xx, l.(xx), label=label)
            plot(x, zeros(n), "ko", ms=4)
        else
            semilogy(xx, abs.(l.(xx)), label=label)
        end
    end
    xlabel(L"x")
    ylabel(L"\ell(x)")
    legend(loc="best", frameon=false)
    display(gcf())
end

function interpolant()
    n = 5
    # n = 15
    f = x->1/(1 + 4*x^2)

    clf()
    xx = LinRange(-1,1,1000)
    plot(xx, f.(xx), "k-")
    for (label,x) = (
        ( "uniform", LinRange(-1,1,2n+1)[2:2:end-1] ),
        # ( "Chebyshev", cos.(LinRange(0,π,2n+1)[2:2:end-1]) ),
    )
        plot(xx, interpolate.(Ref(x), Ref(f.(x)), xx), label=label)
        plot(x, f.(x), "ko", ms=4)
    end
    xlabel(L"x")
    ylabel(L"\ell(x)")
    display(gcf())
end


# Helper functions

chebyshev_points(n) = cos.(LinRange(0,π,2n+1)[2:2:end-1])
function sample_chebyshev_interpolant(f,n,k=10)
    # Input:
    #    f : function to interpolate
    #    n : number of interpolation points
    #    k : Oversampling rate
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

    nn = k*n
    xx = chebyshev_points(nn)
    pxx = FFTW.r2r([c; zeros(nn-n)], FFTW.REDFT01)

    return xx,pxx
end

function chebyshev_interpolation_error(f,n)
    xx,pxx = sample_chebyshev_interpolant(f,n)
    return norm(pxx .- f.(xx), Inf)
end

function convergence_abs()
    n = 2 .^ (1:10)
    f = abs
    errors = chebyshev_interpolation_error.(f,n)

    clf()
    loglog(n, errors, "-o", ms=2)
    nn = (8,n[end])
    loglog(nn, 1.2.*nn.^-1,"k--")
    xlabel(L"n")
    ylabel(L"\|f - p\|_{[-1,1]}")
    display(gcf())
end

function convergence_abssin3()
    n = 2 .^ (1:10)
    f = x -> abs(sin(4π*x))^3
    errors = chebyshev_interpolation_error.(f,n)

    clf()
    loglog(n, errors, "-o", ms=2)
    nn = (2^4.5,n[end])
    loglog(nn, 1e4.*nn.^-3,"k--")
    xlabel(L"n")
    ylabel(L"\|f - p\|_{[-1,1]}")
    display(gcf())
end

function plot_abssin3()
    n = 16
    # n = 32
    # n = 64
    x = chebyshev_points(n)
    f = x -> abs(sin(4π*x))^3
    xx,pxx = sample_chebyshev_interpolant(f,n, 2^10÷n)

    clf()
    plot(xx, f.(xx), "k", lw=0.5, label=L"f(x)")
    plot(xx, pxx, label=L"p(x)")
    plot(x, f.(x), "ko", ms=3)
    xlabel(L"x")
    legend(loc="center left", frameon=false, bbox_to_anchor=(1,0.5))
    display(gcf())
end

function convergence_expinvx()
    n = 2 .^ (1:10)
    f = x -> exp(-inv(abs(x)))
    errors = chebyshev_interpolation_error.(f,n)

    clf()
    if true
        loglog(n, errors, "-o", ms=2)
    else
        semilogy(n, errors, "-o", ms=2)
    end
    xlabel(L"n")
    ylabel(L"\|f - p\|_{[-1,1]}")
    display(gcf())
end

function convergence_sin()
    n = 1:100
    f = x -> sin(4π*x)
    errors = chebyshev_interpolation_error.(f,n)

    clf()
    semilogy(n, errors, "-o", ms=2)
    xlabel(L"n")
    ylabel(L"\|f - p\|_{[-1,1]}")
    display(gcf())
end

function convergence_runge()
    n = 1:100
    a = (4,25)
    c = ("C0", "C1")
    s = (0.25,8)
    nn = ((1,70), (5,100))

    clf()
    for (a,c,s,nn) in zip(a,c,s,nn)
        f = x -> inv(1 + a*x^2)
            r = BernsteinEllipses.radius(im/sqrt(a))
        errors = chebyshev_interpolation_error.(f,n)
        semilogy(n, errors, "$c-o", ms=2)
        semilogy(nn, s.*r.^.-nn, "$c--")
    end
    xlabel(L"n")
    ylabel(L"\|f - p\|_{[-1,1]}")
    display(gcf())
end