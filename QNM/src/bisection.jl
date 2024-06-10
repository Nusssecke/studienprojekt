using QNM

function real_wrapper(omega)
    q, k, omega = [1.19225, 0.5, 0.0 + omega * im]
    return real(computeShearMode(q, k, omega))
end

function imag_wrapper(omega)
    q, k, omega = [1.19225, 0.5, 0.0 + omega * im]
    return imag(computeShearMode(q, k, omega))
end

# Find QNM modes
function findQNMBisection(func_real::Function, func_imag::Function, a::Number, b::Number, tol::AbstractFloat=1e-5, maxiter::Integer=100)
    # Compute the function value on both boundaries
    f_a, f_b = [func_real(a), func_real(b)]

    # Check if there is a sign flip in the interval
    if f_a * f_b <= 0 && false
        # There is a sign flip, do bisection on normal interval
        return bisection(func_real, a, b, tol=tol, maxiter=maxiter)
    else
        println("Possibly more than one root in interval!")
        # Sample the function on the interval to find possible roots
        numSamples = 2
        step = (b-a) / numSamples
        samples = [func_real(x) for x in (a+step):step:(b-step)]

        pushfirst!(push!(samples, f_b), f_a) # Add start and end to the array
        println("Sampling; #: $(numSamples-1), step: $step, samples: $samples")
        
        # For each pair check if they have a sign change
        roots = Float64[] # Initialize an empty array to store the results

        # Iterate through the list and check for sign changes
        for i in 1:length(samples)-1
            if samples[i] * samples[i+1] < 0
                root = bisection(func_real, a + step*(i-1), a + step*i, tol=tol, maxiter=maxiter)
                push!(roots, root)
            end
        end
        println("Found roots: $roots")
        println("Checking against imaginary part")

        qnm = []
        # Check if found roots are also roots for the imaginary part
        epsilon = 0.001
        for root in roots
            imag_root = bisection(func_imag, root - epsilon, root + epsilon, tol=tol, maxiter=maxiter)
            println("Imag root $imag_root")
            if !isempty(imag_root)
                push!(qnm, root)
            end
        end
        return qnm
    end
end

# Serarch for roots in the real part and check in the imaginary part if true root
function bisection(f::Function, a::Number, b::Number; tol::AbstractFloat=1e-5, maxiter::Integer=100)
    f_a, f_b = [f(a), f(b)]
    f_a * f_b <= 0 || error("No real root in [a, b]") # Error because we checked before
    
    iterations = 0
    local c
    while (b-a) > tol
        iterations += 1
        iterations != maxiter || error("Max iteration exceeded")
        c = (a+b)/2 # Midpoint
        f_c = f(c)

        if f_c == 0
            break
        elseif f_a * f_c > 0
            a = c  # Root is in the right half of [a,b].
            f_a = f_c
        else
            b = c  # Root is in the left half of [a,b].
        end
    end

    return c
end

# There is a root inbetween -0.7...-0.65
qnm = findQNMBisection(real_wrapper, imag_wrapper, -5.0, -0.5)
println(qnm)


# function aFunction(x::Float64)
    # return x^4-x^2
# end
# roots = findQNMBisection(aFunction, aFunction, 4.0, 7.0)
# println(roots)