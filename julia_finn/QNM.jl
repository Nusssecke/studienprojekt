using DifferentialEquations
using BenchmarkTools
using BoundaryValueDiffEq
using Plots
using NLsolve

include("horizon_expansion.jl")

f(u, qt) = 1 - (1+qt^2) * u^2 + qt^2 * u^3 # blackening_factor
df(u, qt) = -2 * (1+qt^2) * u + 3 * qt^2 * u^2 # derivative of blackening_factor


function shear_mode_eq!(du, u, p, t)
    # u is current state variable, du is the derivative of u at time t, t is current time
    c0, qt, omega, kk = p # p is a vector of parameters

    phiReal, phiImag, dphiReal, dphiImag = u
    du[1] = dphiReal
    du[2] = dphiImag
    phi = phiReal + im * phiImag
    dphi = dphiReal + im * dphiImag
    eq = (
            -((omega^2 - kk^2*(-1 + t)*(-1 + t*(-1 + qt^2*t)))* phi)/(4*(-1 + t)^2*t*(1 + t - qt^2*t^2)^2) 
            -((-1 + t*(-1 + qt^2*t))*(-1 + t^2*(-1 + qt^2*(-1 + 2*t)))* dphi)/((-1 + t)*t*(1 + t - qt^2*t^2)^2)
        ) * (1 + t - qt^2*t^2)^2/((-1 + t*(-1 + qt^2*t))^2)
    du[3] = real(eq)
    du[4] = imag(eq)
    nothing
end

# boundary conditions
function boundary_condition!(residual, u, p, t)
    residual[1] = real(phiHorizonExpansion(uHorizonNumerical, p...)) - u[end][1] # solution at horizon should be the horizon expansion
    residual[2] = imag(phiHorizonExpansion(uHorizonNumerical, p...)) - u[end][2] # solution at horizon should be the horizon expansion
    residual[3] = real(dphiHorizonExpansion(uHorizonNumerical, p...)) - u[end][3] # derivative at horizon should be the derivative of the horizon expansion
    residual[4] = imag(dphiHorizonExpansion(uHorizonNumerical, p...)) - u[end][4] # derivative at horizon should be the derivative of the horizon expansion
end

epsilon = big(1.0)/big(1000000000.0)
uBoundaryNumerical = epsilon
uHorizonNumerical = big(1.0)-big(1.0)/big(10.0)

function shear_mode(qt, omega, kk, c0=1)
    # println("Starting calculation of phiSol with qt = ", qt, ", omega = ", omega, ", kk = ", kk)

    u0 = [1, 0, 1, 0]; # initial conditions    
    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
    parameters = [c0, qt, omega, kk]

    boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)
    # dtmax = 1/1000, reltol = 1e-12, abstol = 1e-12
    solution = solve(boundary_value_problem, MIRK2(), dt= 0.05)

    # Plot the first to variables of the solution atop one another
    display(plot(solution))
    # println("Finished calculation of phiSol with qt = ", qt, ", omega = ", omega, ", kk = ", kk)

    return u -> [solution(u, idxs = 1), solution(u, idxs = 2)]
end

# sol = phiSol(0, 1/2, 1/2)
# println(sol(uBoundaryNumerical))
# println(sol(1/2))
sol = shear_mode(sqrt(2) - 1/100, 0.0 + im * -0.003, 1)

function solve_quasinormal_mode(omega)
    println("Omega = ", omega)
    solution =  shear_mode(sqrt(2) - 1/100, omega[1] + im * omega[2], 1)
    println("Solution = ", solution(uBoundaryNumerical))
    return solution(uBoundaryNumerical)
end

# qnmRoutine
function qnmRoutine(initial_guess_omega)
    s = nlsolve(solve_quasinormal_mode, initial_guess_omega, xtol=1e-30, ftol=1e-30)
    println(s.zero)
end

# qnmRoutine([0.0, -0.01])

#=
values = zeros(Complex{Float64}, 10, 10)
for i in 1:10
    for j in 1:10
        println("Omega = ", i/100 - im * j/100)
        value = (shear_mode(sqrt(2) - 1/100 , i/100 - im * j/100, 1))(uBoundaryNumerical)
        println("Value = ", value[1], " + im * ", value[2])
        values[i, j] = value[1] + im * value[2]
    end
end

print(values)
=#