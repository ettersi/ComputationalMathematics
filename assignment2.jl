using PyPlot


##########
# Task 2 #
##########

"""
hermite_interpolate(f,x)

Evaluate `p(x)` where `p` is the unique polynomial specified by the
interpolation problem on the assignment sheet.
"""
function hermite_interpolate(f,x)
    # TODO: Your code here
    return NaN
end

function draw_heart()
    f = [
          0.0   0.0
          0.2   2.0

          1.0   0.0
          0.0  -1.8

          0.0  -1.2
         -1.0  -1.5
    ]

    t = LinRange(0,1,1000)
    clf()
    for s in (+1,-1)
        plot(
            s .* hermite_interpolate.(Ref(f[1:4,1]),t),
                 hermite_interpolate.(Ref(f[1:4,2]),t),
        )
        plot(
            s .* hermite_interpolate.(Ref(f[3:6,1]),t),
                 hermite_interpolate.(Ref(f[3:6,2]),t),
        )
    end
    axis("equal")
    display(gcf())
end



##########
# Task 3 #
##########

using FastGaussQuadrature
function composite_gauss(f,a,b,m,n)
    # TODO: Your code here
    return NaN
end

function composite_gauss_convergence()
    f = sin
    a,b = 0,π
    Iref = 2

    clf()
    for (i,n) = enumerate(1:3)
        d = 2n
        m = round.(Int, 10.0.^LinRange(0.5,3,20))
        errors = abs.(composite_gauss.(f,a,b,m,n) .- Iref)
        loglog(m, errors, "C$(i-1)", label=latexstring("n = $n"))
        s = 2*errors[1]*m[1]^d
        loglog(m, s.*float.(m).^(-d), "C$(i-1)--", label=latexstring("O(m^{-$d})"))
    end
    ylim([1e-16,1])
    legend(frameon=false, loc="center left", bbox_to_anchor=(1,0.5))
    xlabel(L"m")
    ylabel(L"\mathrm{error}(m)")
    display(gcf())
end
