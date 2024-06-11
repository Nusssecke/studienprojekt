using SymEngine, QuasinormalModes, RootsAndPoles, Plots, LsqFit, CSV, DataFrames, Base.Threads

# ------------------------------------------------------------------
# Numeric Reissner-Nordström Black Hole
# ------------------------------------------------------------------

struct NReissnerNordströmData{N, T} <: NumericAIMProblem{N, T}
	nIter::N
	x0::T
	q::Float64 # charge
	k::Float64 # wave number
end

function NReissnerNordströmData(nIter::N, x0::T, q::Float64, k::Float64) where {N, T}
	return NReissnerNordströmData{N, T}(nIter, x0, q, k)
end

QuasinormalModes.λ0(d::NReissnerNordströmData{N, T}) where {N, T} =
	(u, w) -> ((#((1/((-1+u)^2*u)/(4*(1 + u - d.q^2*u^2)^2)))*(- (4*(-1 + u)*(-1 + u*(-1 + d.q^2*u))*(-1 + u^2*(-1 + d.q^2*(-1 + 2*u))) +  
				#(4*(-1 + u)*(1 + u - d.q^2*u^2)^2*(2*(-2 + d.q^2)*(-1 + u) + im*u*w))/(-2 + d.q^2))))
		-(4*(-1 + u)*(-1 +u*(-1 + d.q^2*u))*((-2 + d.q^2)*(1 +  u^2*(-3 + d.q^2*(-3 + 4*u))) + 
	    im*u*(-1 + u*(-1 + d.q^2*u))*w))/(-2 + d.q^2))/(4*(-1 + u)^2*u*(1 + u - d.q^2*u^2)^2))
QuasinormalModes.S0(d::NReissnerNordströmData{N, T}) where {N, T} =
	(u, w) -> (#((1/((-1+u)^2*u)/(4*(1 + u - d.q^2*u^2)^2)))*(- (-d.k^2)*(-1 + u)*(-1 + u*(-1 + d.q^2*u)) + 	w^2 + 
				#(2*(-1 + u*(-1 + d.q^2*u))*(-1 + u^2*(-1 + d.q^2*(-1 + 2*u)))*   (2*(-2 + d.q^2)*(-1 + u) + im*u*w))/((-2 + d.q^2)*u) - ((1 + u - d.q^2*u^2)^2*w*(-2*
		  		#im*(-2 + d.q^2)*(-2 + u) + u*w))/(-2 + d.q^2)^2))
			(-(((-(((-1 + u)*(-1 + u*(-1 + d.q^2*u))*(4 + u*(d.k^2 + 4*u*(1 + d.q^2 - 2*d.q^2*u))))/u) + w^2 + (2 *im*(-1 + u)*(-1 + u*(-1 + d.q^2*u))*(-1 + 
	  u*(-2 + 3*d.q^2*u))*w)/(-2 + d.q^2) - (u*(1 + u - d.q^2*u^2)^2*w^2)/(-2 + d.q^2)^2))
	)/(4*(-1 + u)^2*u*(1 + u - d.q^2*u^2)^2)
	))

QuasinormalModes.get_niter(d::NReissnerNordströmData{N, T}) where {N, T} = d.nIter
QuasinormalModes.get_x0(d::NReissnerNordströmData{N, T}) where {N, T} = d.x0

# ------------------------------------------------------------------
# 2. Constructing problems and caches
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# 3. RootsAndPoles.jl search domain and mesh construction
# ------------------------------------------------------------------
params = GRPFParams(
	100,     # the maximum number of refinement iterations before `grpf` returns.
    50000000,   # the maximum number of Delaunay tessalation nodes before `grpf` returns.
    5,       # maximum ratio of the longest to shortest side length of Delaunay triangles before they are split during `grpf` refinement iterations.
    5000,    # provide a size hint to the total number of expected nodes in the Delaunay tesselation. Setting this number approximately correct can improve performance
    1.0e-9, # maximum allowed edge length of the tesselation defined in the `origcoords` domain before returning
	true    # use `Threads.@threads` to run the user-provided function `fcn` across the `DelaunayTriangulation`
)
# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------
# Initialize arrays
k_values = 0.0:0.1:3.0
omega_values = []

# Parallele Schleife
Threads.@threads for k in k_values
    m = Serial()
    p = NReissnerNordströmData(UInt(48), Complex(0.5, 0.0),0.99*sqrt(2), k)
    c = AIMCache(p)

    # Define δ_wrapper function
    function δ_wrapper(z)
        computeDelta!(m, p, c, z)
    end

    # Define search domain and mesh
    origcoords = rectangulardomain(Complex(0.0, -5.0), Complex(1.0, -0.0), 1e-3)

    # Define search settings
    params = GRPFParams()

    # Calculate roots
    roots, poles = grpf(δ_wrapper, origcoords, params)
    push!(omega_values, roots)
     
# Extract roots and flatten omega_values
	flat_omega_values = vcat([roots for roots in omega_values]...)

# Create corresponding k_values
	flat_k_values = vcat([fill(k, length(roots)) for (k, roots) in zip(k_values, omega_values)]...)

# Convert the eigenvalues and k values to a DataFrame
	df = DataFrame(K_values = flat_k_values, Eigenvalues = flat_omega_values)

# Write the DataFrame to a CSV file
	CSV.write("99*sqrt2_k1-3.csv", df)
end
#df = CSV.read("99q_k001.csv", DataFrame)

# Extract the k and eigenvalues values
#flat_k_values = df.K_values

# Convert the eigenvalues to complex numbers
#flat_omega_values = parse.(Complex{Float64}, df.Eigenvalues)

#ext_model(x,p)= p[1] .*x.^ p[2] .+ p[3]
# Fit the imaginary part of the eigenmodes to a curve
#result_imag_ext = curve_fit(ext_model, flat_k_values, imag.(flat_omega_values), [1.0,1.0, -3.5])
#println("Fitted parameters for the imaginary part: ", result_imag_ext.param)
#plot(flat_k_values, ext_model(flat_k_values, result_imag_ext.param), label="extended fit for Im(ω)", legend=:outertopleft)
#plot!(flat_k_values ,imag.(flat_omega_values), seriestype = :scatter, xlabel="k", ylabel="Im(ω)",label= "(Im(ω(k)))" ,title="Dispersion \$ω(k) q=0.99q_{extr}\$ ", legend=:outertopleft)#q\$_{extr}\$
println(flat_k_values)
println(imag.(flat_omega_values))
