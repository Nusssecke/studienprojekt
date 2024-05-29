using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

#---------------------------------------------------------#
# Plot QNM #

function plotQNM(omega_real, omega_imag, q, k, plot_graph=false, dtmax=0.1)
    parameters = [1.0, q, k, omega_real + im * omega_imag] # Constants for the differential equation
    expansionValues = [
        real(phiHorizonExpansion14(uHorizonNumerical, parameters...)),
        imag(phiHorizonExpansion14(uHorizonNumerical, parameters...)),
        real(dphiHorizonExpansion14(uHorizonNumerical, parameters...)),
        imag(dphiHorizonExpansion14(uHorizonNumerical, parameters...))
    ]
    parameters = vcat(parameters, expansionValues)
    println("Horizon calculation done")

    # I am still not sure what the real boundary conditions are
    u0 = [0, 0, 0, 0]; # initial conditions
    # u0 = [real(phiHorizonExpansion(uHorizonNumerical, parameters...)), imag(phiHorizonExpansion(uHorizonNumerical, parameters...)), real(dphiHorizonExpansion(uHorizonNumerical, parameters...)), imag(dphiHorizonExpansion(uHorizonNumerical, parameters...))]; # initial conditions
	tspan = (QNM.uBoundaryNumerical, QNM.uHorizonNumerical) # Possible values for u
    boundary_value_problem = BVProblem(simple_shear_mode_eq!, simple_boundary_condition!, u0, tspan, parameters)

    solution = solve(boundary_value_problem, Shooting(Rodas5P()), dt=dtmax/10, dtmax=dtmax)

    if plot_graph
        display(plot(solution, title="QNM"))
        # display(plot(solution, plotdensity = 10000, idxs=[(1,2,0), (3,4,0)], xaxis="R(phi)", yaxis="I(phi)"))
        # display(plot(solution, plotdensity = 10000, idxs=[(1,2),(3,4)]))
        # https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/
    end

end
#---------------------------------------------------------#

println("Plot Test")
plotQNM(3.0, -3.0, sqrt(2)*0/((3)^(3/4)), 0.0, true)