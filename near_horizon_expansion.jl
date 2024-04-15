using Symbolics
@variables w q alpha u

# Find near horizon expansion

alpha = im * w / (2 * (-2 + q^2))
println(substitute(alpha, q => 1))

f = (1 - u)

# Function for computing the nth coefficient
function coefficient_horizon(n)
    println("Computing coefficient $n")

end