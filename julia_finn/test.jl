include("horizon_expansion.jl")

print("Test Case 1\n")
u = 1/2
c0 = 1
qt = 1/10
omega = 1/10
kk = 0
println(phiHorizonExpansion(u, c0, qt, omega, kk))
# Expected Output: 0.998688687046780895747482 + 0.0246516649076105618329428 i

print("Test Case 2\n")
println(dphiHorizonExpansion(u, c0, qt, omega, kk))
# Expected Output: 0.00123882900990351218987474 + 0.0344982622722421530720753 I

using NOMAD

function eval_fct(omegas)
    println("eval_fct with omegas = ", omegas[1], " ", omegas[2], "im")
    omega = omegas[1] + im * omegas[2]
    a_sol = phiSol(0, omega, 1)
    output = a_sol(uBoundaryNumerical)
    println("output = ", output)
    # constraint imaginary part smaller than zero
    bb_outputs = [output[1]^4 + output[2]^4; output[2]]
    success = true
    count_eval = true
    return (success, count_eval, bb_outputs)
end

pb = NomadProblem(2, # number of inputs of the blackbox
                  2, # number of outputs of the blackbox
                  ["OBJ", "EB"], # type of outputs of the blackbox
                  eval_fct;
                  lower_bound=[0, -5.0],
                  upper_bound=[5.0, 0])

result = NOMAD.solve(pb, [4.0, -3.0])