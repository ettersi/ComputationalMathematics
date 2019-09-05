using CSV
using LinearAlgebra

# Load the data
# Dataset available from https://vincentarelbundock.github.io/Rdatasets/csv/ggplot2/diamonds.csv
data = CSV.read("diamonds.csv");
n = size(data,1)
clarity_codes = ("I2","I1","SI2","SI1","VS2","VS1","VVS2","VVS1","IF")
colour_codes = ("J","I","H","G","F","E","D")
cut_codes = ("Fair","Good","Very Good","Premium","Ideal")
clarity_offset = 2
colour_offset = clarity_offset + length(clarity_codes)
cut_offset = colour_offset + length(colour_codes)

# Assemble matrix and rhs
A = hcat( # Horizontally concatenate vectors passed as arguments
    ones(n),
    Float64.(data[!,:carat]),
    ((c->data[!,:clarity] .== c).(clarity_codes))...,
    ((c->data[!,:color] .== c).(colour_codes))...,
    ((c->data[!,:cut] .== c).(cut_codes))...
)
b = Float64.(data[!,:price])

# Assemble indices for all columns except the ones corresponding to c_star and p_star
clarity_idx = findfirst(isequal("SI1"), clarity_codes)
colour_idx = findfirst(isequal("H"), (colour_codes))
cut_idx = findfirst(isequal("Ideal"), (cut_codes))
idx = [
    1;2;
    clarity_offset .+ setdiff(1:length(clarity_codes),clarity_idx);
    colour_offset  .+ setdiff(1:length(colour_codes),colour_idx)
    cut_offset     .+ setdiff(1:length(cut_codes),cut_idx)
]

# Solve the LSQ problem
x = zeros(size(A,2))
x[idx] = A[:,idx]\b

# Print the results
println("Offset: ", round(x[2],digits=2))
println("Price per carat: ", round(x[2],digits=2))
println()
println("Clarity surcharge:")
for i = 1:length(clarity_codes)
    if i < 3
        spaces = "   "
    elseif i < 7
        spaces = "  "
    elseif i < 9
        spaces = " "
    else
        spaces = "   "
    end
    println(spaces,clarity_codes[i],": ",round(x[clarity_offset+i],digits=2))
end
println()
println("Colour surcharge:")
for i = 1:length(colour_codes)
    println(" ",colour_codes[i],": ",round(x[colour_offset+i],digits=2))
end
println()
println("Cut surcharge:")
for i = 1:length(cut_codes)
    println(" ",cut_codes[i],": ",round(x[cut_offset+i],digits=2))
end
println()
println("Price uncertainty: $(norm(A*x-b)/sqrt(n))")
