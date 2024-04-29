using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

function rescale_to_piT(omega, qt)
    return omega * 2 / (2 - qt^2)
end

function rescale_to_massdensitiy(omega, qt)
    return omega / sqrt( sqrt((1 + qt^2) / 1) )
end

# Initial conditions and parameters
qt_extremal = sqrt(2)

function solveIt(solver)
    omega_real, omega_imag = (big(312)/big(100), -big(275)/big(100))
    c0 = 1
    qt = 1/10
    kk = 0
    u0 = [big(1.0)/big(2.0), big(1.0)/big(2.0), 0, 0, omega_real, omega_imag]; # initial conditions

    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
    parameters = [c0, qt, kk]

    boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)

    # Shooting(Rodas4P()) not good enough
    solution = solve(boundary_value_problem, Shooting(solver), dtmax=0.01)
    println("omega = ", solution.u[1][5] + solution.u[1][6] * im)
    return nothing 
end

# solveIt(Rosenbrock23()) # Quite good 2 Nachkommastellen
# solveIt(Rodas5()) # schnell
solveIt(Rodas5P())
# solveIt(ESDIRK547L2SA2()) ganz schlect