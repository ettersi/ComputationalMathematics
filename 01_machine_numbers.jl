using BenchmarkTools

function default_integer_type()
    @show Int
    @show typeof(42)
end

function binary_representation()
    @show bitstring(Int8(0))
    @show bitstring(Int8(1))
    @show bitstring(Int8(2))
    @show bitstring(Int8(3))
    println()
    @show bitstring(Int8(-1))
    @show bitstring(Int8(-2))
    @show bitstring(Int8(-3))
end

function type_limits()
    @show typemax(Int8)
    @show typemin(Int8)
    println()
    @show typemax(UInt8)
    @show typemin(UInt8)
end

function overflow()
    @show Int8(127) + Int8(1)  # overflow
    @show Int8(-128) - Int8(1) # underflow
end

function conversion_of_unrepresentabe_values()
    try
        Int(0.5)
    catch e
        Base.display_error(e)
    end
    try
        Int8(128)
    catch e
        Base.display_error(e)
    end
    try
        UInt8(-1)
    catch e
        Base.display_error(e)
    end
end

function bigint_performance()
    # Utility function to time addition in type T
    function benchmark(T)
        # The following line is in principle the same as
        #   @btime T(0) + T(0)
        # The extra complication is necessary to prevent the compiler
        # from replacing `0 + 0` -> `0`, which would render the
        # timings meaningless.
        @btime $(Ref(T(0)))[] + $(Ref(T(0)))[]
    end

    println("Int addition:")
    benchmark(Int)
    println()
    println("BigInt addition:")
    benchmark(BigInt)
end

function type_promotion()
    @show typeof(Int8(0) + Int8(0))
    @show typeof(Int8(0) + Int16(0))
    @show typeof(Int8(0) + 0) # Recall `typeof(0) == Int`
end



function div_performance()
    n = 1000
    x = zeros(n)
    y = zeros(n)
    a = 1.0

    @btime $y .= $x./$a       # Divide n times
    @btime $y .= $x.*inv($a)  # Divide once
    # The dots here mean to apply the operations elementwise.
    # Example: `[1,2,3].^2 == [1,4,9]`
    # The $ are required for `@btime`. I recommend to simply ignore them, but
    # if you want to know more then have a look at
    # https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md

    # Division is more expensive than multiplication, so the second version
    # runs faster. The compiler is not allowed to make this optimisation for us
    # because it could change the result.
end