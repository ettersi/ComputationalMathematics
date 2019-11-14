function integers()
    println(bitstring(0)[end-4:end]) # -> "00000"
    println(bitstring(1)[end-4:end]) # -> "00001"
    println(bitstring(2)[end-4:end]) # -> "00010"
    println(bitstring(3)[end-4:end]) # -> "00011"
    println(bitstring(4)[end-4:end]) # -> "00100"
    println(bitstring(5)[end-4:end]) # -> "00101"
    println(bitstring(6)[end-4:end]) # -> "00110"
    println(bitstring(7)[end-4:end]) # -> "00111"
end

function integers_overflow()
    println(typemax(Int))   # ->  9223372036854775807
    println(typemax(Int)+1) # -> -9223372036854775808
end

function floats_sign()
    println(bitstring( 1.0)[1]) # -> "0"
    println(bitstring(-1.0)[1]) # -> "1"
end

function floats_exponent()
    println(bitstring(0.5)[2:12]) # -> "01111111110"
    println(bitstring(1.0)[2:12]) # -> "01111111111"
    println(bitstring(2.0)[2:12]) # -> "10000000000"
    println(bitstring(4.0)[2:12]) # -> "10000000001"
end

function floats_mantissa()
    println(bitstring(1.0 )[13:15]) # -> "000"
    println(bitstring(1.25)[13:15]) # -> "010"
    println(bitstring(1.5 )[13:15]) # -> "100"
    println(bitstring(1.75)[13:15]) # -> "110"
end

function floats_special_values()
    println(bitstring(Inf)[ 2:12]) # -> "11111111111"
    println(bitstring(Inf)[13:15]) # -> "000"
    println(bitstring(NaN)[ 2:12]) # -> "11111111111"
    println(bitstring(NaN)[13:15]) # -> "100"
end

function floats_addition()
    println( (-1.0 +  1.0) + eps()/2  ) # -> 1.1102230246251565e-16
    println(  -1.0 + (1.0  + eps()/2) ) # -> 0.0
end
