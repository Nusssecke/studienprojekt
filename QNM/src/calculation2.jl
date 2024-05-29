using QNM, RootsAndPoles, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations

# ------------------------------------------------------------------
# Numeric Reissner-Nordström Black Hole
# ------------------------------------------------------------------

# TODO TODO TODO Change to the simple boundary value problem
function computeQNM(q, k, omega_real, omega_imag, solver=Rodas5P(), dtmax=0.01)
    parameters = [1.0, q, k, omega_real + im * omega_imag] # Constants for the differential equation
    
    # Pre calculate the values for the boundary conditions
    horizonExpansionValues = [real(phiHorizonExpansion14(uHorizonNumerical, parameters...)), imag(phiHorizonExpansion14(uHorizonNumerical, parameters...)), real(dphiHorizonExpansion14(uHorizonNumerical, parameters...)), imag(dphiHorizonExpansion14(uHorizonNumerical, parameters...))]
    parameters = vcat(parameters, horizonExpansionValues)

    # I am still not sure what the real boundary conditions are
    u0 = [0, 0, 0, 0]; # initial conditions
    # u0 = [real(phiHorizonExpansion(uHorizonNumerical, parameters...)), imag(phiHorizonExpansion(uHorizonNumerical, parameters...)), real(dphiHorizonExpansion(uHorizonNumerical, parameters...)), imag(dphiHorizonExpansion(uHorizonNumerical, parameters...))]; # initial conditions
	tspan = (QNM.uBoundaryNumerical, QNM.uHorizonNumerical) # Possible values for u
    boundary_value_problem = BVProblem(simple_shear_mode_eq!, simple_boundary_condition!, u0, tspan, parameters)

    solution = solve(boundary_value_problem, Shooting(Rodas5P()), dt=dtmax/10, dtmax=dtmax)
    return solution[1][1] + im * solution[2][1]
end

println("Starting sample computation")
q, k, omega_real, omega_imag = [sqrt(2)*0.1/((3)^(3/4)), 0.0, 3.0, -3.0]
omega = computeQNM(q, k, omega_real, omega_imag, Rodas5P(), 0.01)
println("Ending sample computation")
println(omega)

# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------

function δ_wrapper(z)
    q, k, omega_real, omega_imag = [0.0, 0.0, real(z), imag(z)]
    computeQNM(q, k, omega_real, omega_imag)
end

# ------------------------------------------------------------------
# Search domain and mesh construction
# ------------------------------------------------------------------
const xb = 0.0 # real part begin
const xe = 10.0   # real part end
const yb = -10.0 # imag part begin
const ye = -0.0 # imag part end
const r = 1e-3 # initial mesh step
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
    false    # use `Threads.@threads` to run the user-provided function `fcn` across the `DelaunayTriangulation`
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