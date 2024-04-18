using DifferentialEquations
using BoundaryValueDiffEq
using Plots


const g = 9.81
L = 1.0
tspan = (0.0, pi / 2)
function simplependulum!(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end

function f(u, qt) # blackening_factor
    return 1 - (1+qt^2) u^2 + qt^2 u^3
end

function df(u, qt) # derivative of blackening_factor
    return -2 (1+qt^2) u + 3 qt^2 u^2
end

# u is current state variable
# p is all parameters
# t is current time
# du is the derivative of u at time t
function phi(du, u, p, t)
    phi = u[1]
    dphi = u[2]
    du[1] = dphi
    du[2] = dphi * (f(u, qt) - df(u, qt)) / (u * f(u, qt))
        - phi * (omega^2 - k^2 f(u, qt)) / 4 * rH^2 * u * f(u, qt)^2

end

function bc1!(residual, u, p, a)
    residual[1] = u[end ÷ 2][1] + pi / 2 # the solution at the middle of the time span should be -pi/2
    residual[2] = u[end][1] - pi / 2 # the solution at the end of the time span should be pi/2
end

bvp1 = BVProblem(simplependulum!, bc1!, [pi / 2, pi / 2], tspan)
sol1 = solve(bvp1, MIRK4(), dt = 0.05)
plot(sol1)