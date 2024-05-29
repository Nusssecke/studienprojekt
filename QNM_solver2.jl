using SymEngine
using QuasinormalModes
using RootsAndPoles 
using BenchmarkTools
#myphi = -(f(u)-u*f'(u))*p'(u)/(u*f(u))+p''(u)+ p(u)*(w²-k²*f(u))/(4*f(u)²*u)
#f(u, qt) = 1 - (1+q^2) * u^2 + q^2 * u^3 # blackening_factor
#df(u, qt) = -2 * (1+q^2) * u + 3 * q^2 * u^2 # derivative of blackening_factor

# AdS Boundary located at u = 0
# Horizon located at u = 1

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

const m = Serial()
const p = NReissnerNordströmData(UInt(48), Complex(0.5, 0.0),0.0, 0.0);
const c = AIMCache(p)

# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------

function δ_wrapper(z)
    computeDelta!(m, p, c, z)
end

#search domain and mesh construction
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