using QNM # includes the horizonExpansion and the differential equation
using OrdinaryDiffEq # for the ODE solver required by shooting method
using BoundaryValueDiffEq # for the boundary value problem
using Plots # for plotting

function shear_mode(qt, omega_real, omega_imag, kk, c0=1)
    println("Starting calculation of phiSol with qt = ", qt, ", omega = ", omega_real + im * omega_imag, ", kk = ", kk)

    u0 = [0, 0, 0, 0, omega_real, omega_imag]; # initial conditions    
    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
    parameters = [c0, qt, kk]

    boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)
    # dtmax = 1/1000, reltol = 1e-12, abstol = 1e-12
    solution = solve(boundary_value_problem, Shooting(Rodas4P()), dtmax=0.1)

    # Plot the first to variables of the solution atop one another
    # display(plot(solution))
    println("Finished calculation of phiSol with qt = ", qt, ", kk = ", kk)
    # return solution
    return solution
end

sol = shear_mode(sqrt(2) - 1/100, 0.0, -0.003, 1)
println(sol(uBoundaryNumerical, idxs=1), " ", sol(uBoundaryNumerical, idxs=2))
println(sol(uBoundaryNumerical, idxs=5), " ", sol(uBoundaryNumerical, idxs=6))