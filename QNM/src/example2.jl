using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

function rescale_to_piT(omega, qt)
    return omega * 2 / (2 - qt^2)
end

function rescale_to_massdensitiy(omega, qt)
    return omega / sqrt( sqrt((1 + qt^2) / 1) )
end

# Initial conditions and parameters
qt_extremal = sqrt(2)

function solveIt(solver, omega_real, omega_imag, qt, kk, dtmax=0.01)
    c0 = 1
    # I am still not sure what the real boundary conditions are
    # u0 = [big(0.0), big(0.0), 0, 0, omega_real, omega_imag]; # initial conditions

    vars = [c0, qt, omega_real + im * omega_imag, kk]
    u0 = [real(phiHorizonExpansion14(uHorizonNumerical, vars...)), imag(phiHorizonExpansion14(uHorizonNumerical, vars...)), real(dphiHorizonExpansion14(uHorizonNumerical, vars...)), imag(dphiHorizonExpansion14(uHorizonNumerical, vars...)), omega_real, omega_imag]; # initial conditions
    
    # u0 = [0, 0, 0, 0, omega_real, omega_imag]; # initial conditions

    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
    parameters = [c0, qt, kk]

    boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)

    # Shooting(Rodas4P()) not good enough
    solution = solve(boundary_value_problem, Shooting(solver), dtmax=dtmax)
    println("omega = ", solution.u[1][5] + solution.u[1][6] * im)
    plot(solution, vars=(1, 2, 3, 4))
    return nothing 
end

# solveIt(Rosenbrock23()) # Quite good 2 Nachkommastellen
# solveIt(Rodas5()) # schnell
# solveIt(Rodas5P())
# solveIt(ESDIRK547L2SA2()) ganz schlect

for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    solveIt(Rodas5P(), 3.0, -3.0, sqrt(2)*q/((3)^(3/4)), 0.0)

end

# for dt in [0.01, 0.001, 0.0001, 0.00001]
#     print(dt)
#     solveIt(Rodas5P(), 3.0, -3.0, sqrt(2)*0/((3)^(3/4)), 0.0, dt)
# end