using QNM, NonlinearSolve

function shearMode((omega_real, omega_imag))
    println("Started evaluation of shearMode function with: omega_real = $omega_real, omega_imag = $omega_imag")
    q, k, omega = 0.0, 0.0, Complex(omega_real, omega_imag)
    result = computeShearMode(q, k, omega)
    println("Returning: $(real(result)), $(imag(result))")
    return real(result), imag(result)
end

u0 = [3.0, -2.0]
prob = NonlinearProblem(shearMode, u0)
sol = solve(prob)