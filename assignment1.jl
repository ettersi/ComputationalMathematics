using PyPlot
using LinearAlgebra
using Random

#############
# Question 3

function expm1_errors()
    # TODO: Sample x over an appropriate range
    # Hint: Have a look at Julia's `LinRange` function.
    x = NaN
    f = x -> exp(x) - 1
    relerror = @. abs(f(x) - f(big(x))) / abs(f(big(x)))

    clf()
    # TODO: Plot `(x,relerror)` to show `relerror = O(inv(x))`.
    display(gcf())
end


#############
# Question 4

function linear_system_error()
    # Assemble the linear system
    Random.seed!(42)
    n = 20
    xx = LinRange(-1,1,n)
    A = xx.^(0:n-1)'
    b = rand(n)

    # Compute the solution with `Float64` and `BigFloat` accuracy
    x = x_F64 = A\b
    x_big = big.(A) \ big.(b)

    # Compute estimated and exact relative errors
    C = NaN # TODO: Your code here
    relerror = NaN # TODO: Your code here

    println("Estimated relative error: ", round(C * eps(), sigdigits=3))
    println("    Exact relative error: ", round(relerror, sigdigits=3))
end
