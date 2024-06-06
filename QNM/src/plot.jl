using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

#---------------------------------------------------------#
# Plot QNM #

function plotQNM(omega_real, omega_imag, q, k, plot_graph=false, dtmax=0.1)
    parameters = [1.0, q, k, omega_real + im * omega_imag] # Constants for the differential equation
    horizonExpansionValues = [phiHorizonExpansion(uHorizonNumerical, parameters...), dphiHorizonExpansion(uHorizonNumerical, parameters...)]
    parameters = vcat(parameters, horizonExpansionValues)

    u0 = [0, 0, 0, 0]; # Initial conditions
	tspan = (QNM.uBoundaryNumerical, QNM.uHorizonNumerical) # Possible values for u
    boundary_value_problem = BVProblem(simple_shear_mode_eq!, simple_boundary_condition!, u0, tspan, parameters)

    solution = solve(boundary_value_problem, Shooting(Rodas5P()), dt=dtmax/10, dtmax=dtmax)

    if plot_graph
        #display(plot(solution, title="QNM"))
        display(plot(solution, idxs=[1, 2], title="QNM"))
        # display(plot(solution, plotdensity = 10000, idxs=[(1,2,0), (3,4,0)], xaxis="R(phi)", yaxis="I(phi)"))
        # display(plot(solution, plotdensity = 10000, idxs=[(1,2),(3,4)]))
        # https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/
    end

end
#---------------------------------------------------------#

println("Plot Test")
plotQNM(3.1194516, -2.7466757, 0.0, 0.0, true)