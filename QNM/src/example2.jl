using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots

omega_real, omega_imag = (0, -3)
c0 = 1
qt = 1/10
kk = 0
u0 = [0, 0, 0, 0, omega_real, omega_imag]; # initial conditions    
tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u
parameters = [c0, qt, kk]
boundary_value_problem = BVProblem(shear_mode_eq!, boundary_condition!, u0, tspan, parameters)
solution = solve(boundary_value_problem, Shooting(Rodas4P()), dtmax=0.1)
