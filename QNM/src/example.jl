using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots

# julia --project

u0 = [1.0, 1.0, -3.0] # initial guess for ψ, dψ and E

bvp = BVProblem(schroedinger!, bc!, u0, xspan)
sol = solve(bvp, Shooting(Rodas4P()), dtmax = 0.1)

x = sol.t
ψ = [sol.u[i][1] for i in 1:length(sol.u)]
ψ = ψ / sqrt(sum(ψ .^ 2))
plt = plot(x, ψ, xlabel = "x [fm]", ylabel = "ψ(x)", lw = 2)
display(plt)

println("E = ", sol.u[1][3], " MeV")
