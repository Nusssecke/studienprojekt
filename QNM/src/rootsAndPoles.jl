using QNM, RootsAndPoles, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations, BenchmarkTools, Base.Threads

# ------------------------------------------------------------------
# Numeric Reissner-Nordström Black Hole
# ------------------------------------------------------------------

function computeShearMode(q, k, omega, solver=Rodas5P(), dtmax=0.01)
    println("Computing QNM with q = ", q, ", k = ", k, ", omega = ", omega)
    parameters = [1.0, q, k, omega] # Constants for the differential equation
    
    # Pre calculate the values for the boundary conditions
    horizonExpansionValues = [phiHorizonExpansion(uHorizonNumerical, parameters...), dphiHorizonExpansion(uHorizonNumerical, parameters...)]
    # horizonExpansionValues = [phiHorizonExpansion14(uHorizonNumerical, parameters...), dphiHorizonExpansion14(uHorizonNumerical, parameters...)]
    parameters = vcat(parameters, horizonExpansionValues)
    println("Finished horizonExpansionValues: ", horizonExpansionValues)

    # I am still not sure what the real boundary conditions are
    u0 = [0, 0, 0, 0]; # initial conditions
    # u0 = [real(phiHorizonExpansion(uHorizonNumerical, parameters...)), imag(phiHorizonExpansion(uHorizonNumerical, parameters...)), real(dphiHorizonExpansion(uHorizonNumerical, parameters...)), imag(dphiHorizonExpansion(uHorizonNumerical, parameters...))]; # initial conditions
	tspan = (QNM.uBoundaryNumerical, QNM.uHorizonNumerical) # Possible values for u

    boundary_value_problem = BVProblem(simple_shear_mode_eq!, simple_boundary_condition!, u0, tspan, parameters)

    solution = solve(boundary_value_problem, Shooting(solver), dt=dtmax/10, dtmax=dtmax)
    # display(plot(solution))

    # We are checking for roots at the boundary (so the first element of the solution)
    return solution[1][1] + im * solution[1][2]
end

q, k, omega = [0.0, 1/2, 1/2]
omega = computeShearMode(q, k, omega, Rodas5P(), 0.001)
println("Test Root:", omega, "; Trying to match: ", 0.986296694628728507856650+0.083572478812639569093990im)
# omega = computeQNM(1.0, k, omega, Rodas5P(), 0.001)

# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------

function δ_wrapper(omega)
    q, k, omega = [0.0, 0.0, omega]
    computeShearMode(q, k, omega)
end

# ------------------------------------------------------------------
# Search domain and mesh construction
# ------------------------------------------------------------------
const xb = 3.0 # real part begin
const xe = 4.0   # real part end
const yb = -3.0 # imag part begin
const ye = -2.0 # imag part end
const r = 1.0e-1 # initial mesh step
const origcoords = rectangulardomain(Complex(xb, yb), Complex(xe, ye), r)

# ------------------------------------------------------------------
# 4. RootsAndPoles.jl search settings
# ------------------------------------------------------------------

# For details, see https://github.com/fgasdia/RootsAndPoles.jl
params = GRPFParams(
    100,     # the maximum number of refinement iterations before `grpf` returns.
    5000000,   # the maximum number of Delaunay tessalation nodes before `grpf` returns.
    5,       # maximum ratio of the longest to shortest side length of Delaunay triangles before they are split during `grpf` refinement iterations.
    5000,    # provide a size hint to the total number of expected nodes in the Delaunay tesselation. Setting this number approximately correct can improve performance
    1.0e-9, # maximum allowed edge length of the tesselation defined in the `origcoords` domain before returning
    true    # use `Threads.@threads` to run the user-provided function `fcn` across the `DelaunayTriangulation`
)

# ------------------------------------------------------------------
# 5. Root finding and printing
# ------------------------------------------------------------------

benchmark_result = @benchmark grpf($δ_wrapper, $origcoords, $params)

roots, poles = grpf(δ_wrapper, origcoords, params)

println("Roots:")
for root in roots
    println(root)
end

println("-------------------------------------------")
println("Poles:")

for pole in poles
    println(pole)
end

println("-------------------------------------------")
println("Benchmark result:")
println(benchmark_result)

#plot the δ_wrapper function
# Define the range of real and imaginary parts
# re_range = 2.7:0.01:3.2
# im_range = -3.1:0.01:-2.6

# # Compute the function values on the grid
# δ_values = [δ_wrapper(Complex(re, im)) for re in re_range, im in im_range]

# # Extract the real and imaginary parts
# re_δ_values = real.(δ_values)
# im_δ_values = imag.(δ_values)

#print(δ_wrapper(Complex(3.1194515911651477, -2.7466757443405068)))

# plot(re_range, im_range, re_δ_values, st=:surface, c=:viridis, xlabel="Re(z)", ylabel="Im(z)", zlabel="Re(δ(z))", title="δ(z)", )
# plot(re_range, im_range, re_δ_values, st=:surface, c=:viridis, xlabel="Re(z)", ylabel="Im(z)", zlabel="Re(δ(z))", title="δ(z)", size=(800, 600))
# plot(re_range, im_range, im_δ_values, st=:surface, c=:viridis, xlabel="Re(z)", ylabel="Im(z)", zlabel="Im(δ(z))", title="δ(z)", size=(800, 600))
# plot(re_range, im_range, re_δ_values, st=:surface, c=:viridis, xlabel="Re(z)", ylabel="Im(z)", zlabel="Re(δ(z))", title="δ(z)", size=(800, 600), camera=(30, 30))