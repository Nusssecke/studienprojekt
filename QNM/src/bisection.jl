using QNM

function δ_wrapper(omega)
    q, k, omega = [1.19225, 0.5, 0.0 + omega * im]
    real(computeShearMode(q, k, omega))
end

# Serarch for roots in the real part and check in the imaginary part if true root

function bisection(f::Function, a::Number, b::Number; tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a+b)/2
        fc = f(c)
        if fc == 0
            break
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return c
end

bisection(δ_wrapper, -0.7, -0.65)