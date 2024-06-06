using QNM, NLsolve

function shearMode(omega_real, omega_imag)
    println("Started evaluation of shearMode function with: omega_real = $omega_real, omega_imag = $omega_imag")
    q, k, omega = 0.0, 0.0, Complex(omega_real, omega_imag)
    result = computeShearMode(q, k, omega)
    println("Returning: $(real(result)), $(imag(result))")
    return real(result), imag(result)
end

println("Perfect: ", shearMode(3.1194516, -2.7466757))

sol = nlsolve(shearMode, [3.0, -2.0], ftol=1e-30)
sol.zero