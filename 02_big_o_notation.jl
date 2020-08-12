using BenchmarkTools

function sum_ij(A)
    s = 0.0
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            s += A[i,j]
        end
    end
    return s
end

function sum_ji(A)
    s = 0.0
    for j = 1:size(A,2)
        for i = 1:size(A,1)
            s += A[i,j]
        end
    end
    return s
end

function matrix_sum()
    n = 10_000
    A = zeros(n,n)
    @time sum_ji(A)
    @time sum_ij(A)
end



function matrix_product(A,B)
    @assert size(A,2) == size(B,1)
    C = zeros(size(A,1),size(B,2))
    for i = 1:size(A,1)
        for j = 1:size(B,2)
            for k = 1:size(A,2)
                C[i,j] += A[i,k]*B[k,j]
            end
        end
    end
    return C
end

function matrix_product()
    n = 1_000
    A = zeros(n,n)
    B = zeros(n,n)
    @time matrix_product(A,B)
    @time A*B
end



function sum_if(x,y)
    s = 0.0
    for i = 1:length(x)
        if x[i]
            s += y[i]
        end
    end
    return s
end

function branch_prediction()
    # Based on https://stackoverflow.com/q/11227809

    n = 200_000
    x_rand = rand(Bool,n)
    x_sort = sort(x_rand)
    y = rand(n)

    @btime sum_if($x_rand,$y)
    @btime sum_if($x_sort,$y)
    # `sum_if(x,y)` is much faster if `x` is sorted because branch prediction
    # is easier if `x` is of the form `x = [0,0, ..., 0,0,1,1, ... 1,1]`
end