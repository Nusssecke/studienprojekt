using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

function rescale_to_piT(omega, qt)
    return omega * 2 / (2 - qt^2)
end

function rescale_to_massdensitiy(omega, qt)
    return omega / sqrt( sqrt((1 + qt^2) / 1) )
end

# Initial conditions and parameters
qt_extremal = sqrt(2)

function computeQNM(omega_real, omega_imag, qt, kk, solver, dtmax=0.01, plot=false)
    c0 = 1.0
    vars = [c0, qt, omega_real + im * omega_imag, kk] # For the horizon expansion

    # I am still not sure what the real boundary conditions are
    # u0 = [0, 0, 0, 0, omega_real, omega_imag]; # initial conditions
    u0 = [real(phiHorizonExpansion(uHorizonNumerical, vars...)), imag(phiHorizonExpansion(uHorizonNumerical, vars...)), real(dphiHorizonExpansion(uHorizonNumerical, vars...)), imag(dphiHorizonExpansion(uHorizonNumerical, vars...)), omega_real, omega_imag]; # initial conditions

    # This crashes due to memory usage?
    # u0 = [real(phiHorizonExpansion14(uHorizonNumerical, vars...)), imag(phiHorizonExpansion14(uHorizonNumerical, vars...)), real(dphiHorizonExpansion14(uHorizonNumerical, vars...)), imag(dphiHorizonExpansion14(uHorizonNumerical, vars...)), omega_real, omega_imag]; # initial conditions
    println("u0 = ", u0)
    println("u0 (more precission)", [real(phiHorizonExpansion14(uHorizonNumerical, vars...)), imag(phiHorizonExpansion14(uHorizonNumerical, vars...)), real(dphiHorizonExpansion14(uHorizonNumerical, vars...)), imag(dphiHorizonExpansion14(uHorizonNumerical, vars...)), omega_real, omega_imag])
    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
    parameters = [c0, qt, kk] # Constants for the differential equation

    boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)

    solution = solve(boundary_value_problem, Shooting(solver), dt=dtmax/10, dtmax=dtmax)

    if plot
        plot(solution)
    end

    println("omega = ", solution.u[1][5] + solution.u[1][6] * im)
    return solution.u[1][5] + solution.u[1][6] * im
end
computeQNM(3.0, -3.0, sqrt(2)*0.1/((3)^(3/4)), 0.0, Rodas5P(), 0.00001)
# plotQNM(3.0, -3.0, sqrt(2)*0/((3)^(3/4)), 0.0)

# Compute for all q
# for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#     computeQNM(3.0, -3.0, sqrt(2)*q/((3)^(3/4)), 0.0, Rodas5P(), 0.00001)
# end

# Test different maximal time steps
# for dt in [0.01, 0.001, 0.0001, 0.00001]
#     print(dt)
#     computeQNM(3.0, -3.0, sqrt(2)*0/((3)^(3/4)), 0.0, Rodas5P(), dt)
# end