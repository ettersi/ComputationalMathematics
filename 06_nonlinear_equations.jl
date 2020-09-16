function bisection(f, a,b)
    # Check that `[a,b]` is a bracketing interval
    @assert sign(f(a)) != sign(f(b))

    # Ensure `a < b`
    if b < a; a,b = b,a; end

    # Do the bisection
    while b > nextfloat(a)
        m = (b + a)/2
        if sign(f(a)) != sign(f(m))
            a,b = a,m
        else
            a,b = m,b
        end
    end

    # Return the result
    return b
end

using Test
function test_bisection()
    # Generic test
    @test abs(bisection(sin, 3.0, 4.0) - π) <= eps(Float64(π))

    # Test that `f(a) == 0`, `f(b) == 0` are all handled
    # correctly both for `x == 0` and `x != 0`.
    @test abs(bisection(x->x,  0.0, 1.0)) <= nextfloat(0.0)
    @test abs(bisection(x->x, -1.0, 0.0)) <= nextfloat(0.0)
    @test abs(bisection(x->x-1, 1.0, 2.0) - 1) <= eps(1.0)
    @test abs(bisection(x->x-1, 0.0, 1.0) - 1) <= eps(1.0)

    # Test input handling
    @test_throws Exception bisection(sin, 1.0, 2.0)
    @test bisection(sin, 4.0, 3.0) == bisection(sin, 3.0, 4.0)
end



using Printf

function bisection_convergence()
    for k = 0:12
        @printf("error(k = %2.d) = %.10f", k, 2.0^(-k))
        println()
        if mod(k,3) == 0; println(); end
    end
end

function newton_convergence()
    f = sin
    df = cos
    x = big(1.0)

    for k = 0:6
        @printf("error(k = %d) = %.100f", k, abs(x))
        println()
        x -= f(x) / df(x)
    end
end



using PyPlot

function newton_convergence_slow()
    clf()
    for (i,m) = enumerate(2:4)
        n = 100
        f = x -> x^m
        df = x -> m*x^(m-1)
        x = 1.0

        x_hist = zeros(n)
        for k = 1:n
            x_hist[k] = x
            x -= f(x) / df(x)
        end

        nn = [50,n]
        r = 1 - inv(m)
        s = 2 * abs(x_hist[nn[1]-1]) *  r^-nn[1]
        semilogy(0:n-1, abs.(x_hist), "C$(i-1)", label=latexstring("f(x) = x^$m"))
        semilogy(nn, s.*r.^nn, "C$(i-1)--", label=latexstring("O(($(m-1)/$m)^{k})"))
    end
    xlabel(L"k")
    ylabel(L"|x_k - x^\star|")
    legend(frameon=false)
    display(gcf())
end



using LinearAlgebra
using Test

function square_root(w)
    f = x -> [
        x[1]^2 - x[2]^2 - real(w),
            2*x[1]*x[2] - imag(w)
    ]
    df = x-> [
        2x[1]   -2x[2]
        2x[2]    2x[1]
    ]

    x = [real(w),imag(w)]
    for i = 1:20
        if norm(f(x)) < 10*eps()*norm(x)
            return x[1] + x[2]*im
        end
        x -= df(x) \ f(x)
    end
    error("Newton's method did not converge. Final iterate is x = $(x[1] + x[2]*im).")
end

function test_square_root()
    @testset "square_root" begin
        @test square_root(1.0) == 1.0
        @test square_root(2.0) ≈ sqrt(2.0)
        @test square_root(1.0im) ≈ sqrt(1.0im)
        @test square_root(2.0im) ≈ sqrt(2.0im)
        @test square_root(1.0+1.0im) ≈ sqrt(1.0+1.0im)
        @test_throws Exception square_root(-1.0)
    end
end
