using Random
using RandomNumbers.PCG

#############
# Sum example

function sum_serial(N)
    s = 0.0
    for i = 1:N
        s += sin(i)
    end
    return s
end

function sum_parallel_dumb(N)
    s = 0.0
    Threads.@threads for i = 1:N
        s += sin(i)
    end
    return s
end

function sum_parallel(N)
    s = zeros(4)
    Threads.@threads for i = 1:N
        s[Threads.threadid()] += sin(i)
    end
    return sum(s)
end

function sum_example()
    N = 100_000_000

    println("Serial runtime")
    @time sum_serial(N)
    println()

    println("Parallel runtime")
    @time sum_parallel(N)
end


#####################
# Monte Carlo example

function monte_carlo_serial(N)
    Random.seed!(42)
    s = 0.0
    for i = 1:N
        s += rand()
    end
    return s
end

function monte_carlo_parallel_dumb(N)
    Random.seed!(42)
    s = zeros(Threads.nthreads())
    Threads.@threads for i = 1:N
        s[Threads.threadid()] += rand()
    end
    return sum(s)
end

function monte_carlo_parallel(N)
    s = zeros(Threads.nthreads())
    rng = [MersenneTwister(i) for i = 1:Threads.nthreads()]
    @inbounds Threads.@threads for i = 1:N
        s[Threads.threadid()] += rand(rng[Threads.threadid()])
    end
    return sum(s)
end

function monte_carlo_example()
    N = 100_000_0000

    println("Serial runtime")
    @time monte_carlo_serial(N)
    println()

    println("Parallel runtime")
    @time monte_carlo_parallel(N)
end
