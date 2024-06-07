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
    q, k, omega = [0.336133, 0.0, omega]
    computeShearMode(q, k, omega)
end

# ------------------------------------------------------------------
# Search domain and mesh construction
# ------------------------------------------------------------------
const xb = 3.0 # real part begin
const xe = 3.1   # real part end
const yb = -3.0 # imag part begin
const ye = -2.9 # imag part end
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

#--------------------------------------------------------------
# Extract roots and flatten omega_values
# flat_omega_values = vcat([roots for roots in omega_values]...)
# 
# # Create corresponding k_values
# flat_k_values = vcat([fill(k, length(roots)) for (k, roots) in zip(k_values, omega_values)]...)
# 
# # Convert the eigenvalues and k values to a DataFrame
# df = DataFrame(K_values = flat_k_values, Eigenvalues = flat_omega_values)
# 
# # Write the DataFrame to a CSV file
# CSV.write("99q_k001.csv", df)
# 
# #--------------------------------------------------------------
# # Fitting values
# df = CSV.read("99q_k001.csv", DataFrame)
# 
# # Extract the k and eigenvalues values
# flat_k_values = df.K_values
# 
# # Convert the eigenvalues to complex numbers
# flat_omega_values = parse.(Complex{Float64}, df.Eigenvalues)
# 
# ext_model(x,p)= p[1] .*x.^ p[2] .+ p[3]
# # Fit the imaginary part of the eigenmodes to a curve
# result_imag_ext = curve_fit(ext_model, flat_k_values, imag.(flat_omega_values), [1.0,1.0, -3.5])
# println("Fitted parameters for the imaginary part: ", result_imag_ext.param)
# plot(flat_k_values, ext_model(flat_k_values, result_imag_ext.param), label="extended fit for Im(ω)", legend=:outertopleft)
# plot!(flat_k_values ,imag.(flat_omega_values), seriestype = :scatter, xlabel="k", ylabel="Im(ω)",label= "(Im(ω(k)))" ,title="Dispersion \$ω(k) q=0.99q_{extr}\$ ", legend=:outertopleft)#q\$_{extr}\$
# println(flat_k_values)
# println(imag.(flat_omega_values))
#--------------------------------------------------------------


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