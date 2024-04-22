using DifferentialEquations
using BoundaryValueDiffEq
using Plots

include("horizon_expansion.jl")

f(u, qt) = 1 - (1+qt^2) * u^2 + qt^2 * u^3 # blackening_factor
df(u, qt) = -2 * (1+qt^2) * u + 3 * qt^2 * u^2 # derivative of blackening_factor


function phi(du, u, p, t)
    # u is current state variable, du is the derivative of u at time t, t is current time
    c0, qt, omega, kk = p # p is a vector of parameters

    phiReal, phiImag, dphiReal, dphiImag = u
    du[1] = dphiReal
    du[2] = dphiImag
    phi = phiReal + im * phiImag
    dphi = dphiReal + im * dphiImag
    eq = (-((omega^2 - kk^2*(-1 + t)*(-1 + t*(-1 + qt^2*t)))* phi)/(4*(-1 + t)^2*t*(1 + t - qt^2*t^2)^2) + -((-1 + t*(-1 + qt^2*t))*(-1 + t^2*(-1 + qt^2*(-1 + 2*t)))* dphi)/((-1 + t)*t*(1 + t - qt^2*t^2)^2))* (1 + t - qt^2*t^2)^2/((-1 + t*(-1 + qt^2*t))^2)
    du[3] = real(eq)
    du[4] = imag(eq)
end

# boundary conditions
function bc1!(residual, u, p, t)
    c0, qt, omega, kk = p
    residual[1] = real(phiHorizonExpansion(uHorizonNumerical, c0, qt, omega, kk)) - u[end][1] # solution at horizon should be the horizon expansion
    residual[2] = imag(phiHorizonExpansion(uHorizonNumerical, c0, qt, omega, kk)) - u[end][2] # solution at horizon should be the horizon expansion
    residual[3] = real(dphiHorizonExpansion(uHorizonNumerical, c0, qt, omega, kk)) - u[end][3] # derivative at horizon should be the derivative of the horizon expansion
    residual[4] = imag(dphiHorizonExpansion(uHorizonNumerical, c0, qt, omega, kk)) - u[end][4] # derivative at horizon should be the derivative of the horizon expansion
end

epsilon = 1/1000000000
uBoundaryNumerical = epsilon
uHorizonNumerical = 1-1/10

function phiSol(qt,omega, kk, c0=1)
    println("Starting calculation of phiSol with qt = ", qt, ", omega = ", omega, ", kk = ", kk)
    u0 = [1, 0, 1, 0]; # initial conditions    
    tspan = (uBoundaryNumerical, uHorizonNumerical) # Possible values for u

    # c0, qt, omega, kk
    p = [1, qt, omega, kk] # parameters

    bvp = BVProblem(phi, bc1!, [1/2, 0, 1/2, 0], tspan, p)
    sol = solve(bvp, Shooting(Vern7()), dtmax = 1/1000, reltol = 1e-12, abstol = 1e-12)

    # Plot the first to variables of the solution atop one another
    plt = plot(sol)
    display(plt)
    solution = u -> [sol(u, idxs = 1), sol(u, idxs = 2)]
    println("Finished calculation of phiSol with qt = ", qt, ", omega = ", omega, ", kk = ", kk)
    return solution
end

# sol = phiSol(0, 1/2, 1/2)
# println(sol(uBoundaryNumerical))
# println(sol(1/2))

# qnmRoutine
using NLsolve
function salv(omega)
    println(omega)
    sol =  phiSol(sqrt(2), omega[1] + im * omega[2], 1)
    println("solution=", sol(uBoundaryNumerical))
    return sol(uBoundaryNumerical)
end

# s = nlsolve(salv, [4.0, -3.0], xtol=1e-30, ftol=1e-30)
# println(s.zero)

for i in 1:100
    println((phiSol(sqrt(2), i/10 - im * i/10, 1))(uBoundaryNumerical))
end