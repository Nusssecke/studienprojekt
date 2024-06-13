using QNM, Plots, Base.Threads, OrdinaryDiffEq, BoundaryValueDiffEq, DifferentialEquations

# Find QNM modes
function findQNMBisection(func_real::Function, func_imag::Function, a::Number, b::Number, numSamples::Integer=10, tol::AbstractFloat=1e-5, maxiter::Integer=100)
    # Compute the function value on both boundaries
    f_a, f_b = [func_real(a), func_real(b)]

    # Check if there is a sign flip in the interval
    if f_a * f_b <= 0 && false
        # There is a sign flip, do bisection on normal interval
        return bisection(func_real, a, b, tol=tol, maxiter=maxiter)
    else
        println("Sampling $numSamples points in the interval [$a, $b]")
        # Sample the function on the interval to find possible roots
        step = (b-a) / numSamples

        # samples = [func_real(x) for x in (a+step):step:(b-step)]
        samples = Vector{Float64}(undef, numSamples-1)
        Threads.@threads for x in 1:numSamples-1
            samples[x] = func_real(a + step*x)
        end

        pushfirst!(push!(samples, f_b), f_a) # Add start and end to the array
        println("Sampling; #: $(numSamples-1), step: $step, samples: $samples")
        
        # For each pair check if they have a sign change
        roots = Float64[] # Initialize an empty array to store the results

        # Iterate through the list and check for sign changes
        for i in 1:length(samples)-1
            if samples[i] * samples[i+1] < 0
                println("Samples have different sign: $(a+step*(i-1)) and $(a+step*(i)); $(samples[i]) and $(samples[i+1])")
                root = bisection(func_real, a + step*(i-1), a + step*i, tol=tol, maxiter=maxiter)
                push!(roots, root)
            end
        end
        println("Found roots: $roots")
        println("Checking against imaginary part")

        qnm = AbstractFloat[]
        # Check if found roots are also roots for the imaginary part
        epsilon = 0.001
        for root in roots
            imag_root = bisection(func_imag, root - epsilon, root + epsilon, tol=tol, maxiter=maxiter)
            if !isempty(imag_root)
                push!(qnm, root)
            end
        end
        return qnm
    end
end

# https://mmas.github.io/bisection-method-julia
# Serarch for roots in the real part and check in the imaginary part if true root
function bisection(f::Function, a::Number, b::Number; tol::AbstractFloat=1e-5, maxiter::Integer=100)
    f_a, f_b = [f(a), f(b)]
    if f_a * f_b >= 0
        println("Before error $(f_a * f_b <= 0)")
        error("No root in [a, b], [$a, $b], [$f_a, $f_b]; Product [$(f_a*f_b)]") # Error because we checked before
    end
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

# ------------------------------------------------------------------
# Calculation
# println("Number of threads: $(Threads.nthreads())")

q = 1.19225
samples = 100

# Calculate spectrum
# found_qnm = Dict()
# Threads.@threads for k in 0:0.1:1.0
#     println("Computing QNM for k = $k:")
#     omega = findQNMBisection(omega->real(computeShearMode(q,k,omega*im,Rodas5P(),0.001)), omega->imag(computeShearMode(q,k,omega*im,Rodas5P(),0.001)), -5.0, 0.0, samples)
#     println("Found QNM: $omega")
#     found_qnm[k] = omega
# end
# println(found_qnm)
# QNM.save_results("99", q, found_qnm)

k = 0.005
# findQNMBisection(omega->real(computeShearMode(q,k,omega*im,Rodas5P(),0.0001)), omega->imag(computeShearMode(q,k,omega*im,Rodas5P(),0.0001)), -4.3, -4.2, samples)

q = sqrt(2)
k = 0.05
findQNMBisection(omega->real(computeShearMode(q,k,omega*im,Rodas5P(),0.01,true)), omega->imag(computeShearMode(q,k,omega*im,Rodas5P(),0.01,true)), -4.3, -4.2, samples)

# ------------------------------------------------------------------
# Test functions
# function aFunction(x::Float64)
#     return x^2
# end
# roots = findQNMBisection(aFunction, aFunction, 1.0, 2.0, 40)
# println(roots)