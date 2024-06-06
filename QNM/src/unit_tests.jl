using QNM
c0 = 1

#---------------------------------------------------------#
# DifferentialEquation
# t, phi, dphi, q, k, omega = [0.5, 0.5, 0.5, 1/10, 0, 1/10]
# println("ddphi: 0==", ddphi(t, phi, dphi, q, k, omega) - 1.6649897854242324)
#---------------------------------------------------------#
# Horizon expansion tests
u, q, k, omega = [1/2, 1/10, 0, 1/10]
println("HorizonExpansion: ", phiHorizonExpansion(u, c0, q, k, omega) - (0.998688687046780895747482 + 0.0246516649076105618329428im))
println("HorizonExpansion: ", dphiHorizonExpansion(u, c0, q, k, omega) - (0.00123882900990351218987474 + 0.0344982622722421530720753im))

# u, q, k, omega = [uBoundaryNumerical, 0.0001, 0.0, 0.0 - 10.0im]
# println("HorizonExpansion: ", phiHorizonExpansion(u, c0, q, k, omega))
#---------------------------------------------------------#
# Test for computeShearMode
q, k, omega = [0, 0, 3.1194516 -2.7466757im]
computeShearMode(q, k, omega)