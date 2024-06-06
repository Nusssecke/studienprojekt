using QNM, IntervalArithmetic, IntervalRootFinding, StaticArrays, Plots




# STOPING USAGE PROBABLY WAY TO COMPLICATED!
#https://juliaintervals.github.io/IntervalRootFinding.jl/latest/



# Define the function which this solves.
# For this root finding algorithm we need to return an SVector
function shearMode((omega_real, omega_imag))
    println("Started evaluation of shearMode function with: omega_real = $omega_real, omega_imag = $omega_imag")
    q, k, omega = 0.0, 0.0, Complex(omega_real, omega_imag)
    result = computeShearMode(q, k, omega)
    return SVector(real(result), imag(result))
end

boundaries = Interval(3.0, 4.0) × Interval(-3.0, -2.0)
# rts = roots(shearMode, boundaries)

function test(x)
    println("Started evaluation of test function with: x = $x")
    println("Returning: $((x^2 - 2)^2 * (x^2 - 3))")
    return (x^2 - 2)^2 * (x^2 - 3)
end

function g( (x1, x2, x3) )
    println("Started evaluation of g function with: x1 = $x1, x2 = $x2, x3 = $x3")
    println("Returning: $(x1+x2), $(x2), $(x3)")
    return SVector(x1 + x2,
                   x2,
                   x3
                  )
end

# roots(test, -10..10)
X = -5..5
rts = roots(g, (3..4) × (6..7) × (8..9))