using QNM, RootsAndPoles, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations, BenchmarkTools, Base.Threads

using CSV, DataFrames

# Print number of threads
println("Number of threads: ", Threads.nthreads())

# ------------------------------------------------------------------
# Numeric Reissner-Nordström Black Hole
# ------------------------------------------------------------------

# Definition for compute shear mode

q, k, omega = [0.0, 1/2, 1/2]
omega = computeShearMode(q, k, omega, Rodas5P(), 0.001)
println("Test Root:", omega, "; Trying to match: ", 0.986296694628728507856650+0.083572478812639569093990im)
# omega = computeQNM(1.0, k, omega, Rodas5P(), 0.001)

# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------

function δ_wrapper(omega)
    q, k, omega = [2, 0.05, omega]
    computeShearMode(q, k, omega, Rodas5P(), 0.00001, true)
end

# ------------------------------------------------------------------
# Search domain and mesh construction
# ------------------------------------------------------------------
const xb = 0.0 # real part begin
const xe = 5.0   # real part end
const yb = -5.0 # imag part begin
const ye = 0.0 # imag part end
const r = 1.0e-2 # initial mesh step
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
    1.0e-4, # maximum allowed edge length of the tesselation defined in the `origcoords` domain before returning
    true    # use `Threads.@threads` to run the user-provided function `fcn` across the `DelaunayTriangulation`
)

# ------------------------------------------------------------------
# 5. Root finding and printing
# ------------------------------------------------------------------

roots, poles = grpf(δ_wrapper, origcoords, params)

omega_values = []
push!(omega_values, roots)

println("Roots:")
for root in roots
    println(root)
end

println("-------------------------------------------")
println("Poles:")

for pole in poles
    println(pole)
end

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