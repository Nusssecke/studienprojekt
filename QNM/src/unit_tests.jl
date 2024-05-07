using QNM

u = 1/2
c0 = 1
qt = 1/10
omega = 1/10
kk = 0

println("Test Case 1")
println(phiHorizonExpansion(u, c0, qt, omega, kk))
# Expected Output: 0.998688687046780895747482 + 0.0246516649076105618329428 im

println("Test Case 2")
println(dphiHorizonExpansion(u, c0, qt, omega, kk))
# Expected Output: 0.00123882900990351218987474 + 0.0344982622722421530720753 im
