using QNM
c0 = 1

#---------------------------------------------------------#
# DifferentialEquation
# t, phi, dphi, q, k, omega = [0.5, 0.5, 0.5, 1/10, 0, 1/10]
# println("ddphi: 0==", ddphi(t, phi, dphi, q, k, omega) - 1.6649897854242324)
#---------------------------------------------------------#
# Horizon expansion tests
u, q, k, omega = [1/2, 1/10, 0, 1/10]
println("HorizonExpansion: ", phiHorizonExpansion(q, k, omega) - (0.998688687046780895747482 + 0.0246516649076105618329428im))
println("HorizonExpansion: ", dphiHorizonExpansion(q, k, omega) - (0.00123882900990351218987474 + 0.0344982622722421530720753im))

# u, q, k, omega = [uBoundaryNumerical, 0.0001, 0.0, 0.0 - 10.0im]
# println("HorizonExpansion: ", phiHorizonExpansion(u, c0, q, k, omega))
#---------------------------------------------------------#
# Extremal horizon expansion tests
k, omega = [1/10, 1/10]
println("Extremal HorizonExpansion: ", phiHorizonExpansionExtremal(0, k, omega) - (-0.4594359866445343332230777 + 0.8873686398549546217040529im))
println("Extremal HorizonExpansion: ", dphiHorizonExpansionExtremal(0, k, omega) - (14806.34477565232737106271 + 7666.83821579075083718635im))

#---------------------------------------------------------#
# Test for computeShearMode
q, k, omega = [0, 0, 3.1194516 -2.7466757im]
computeShearMode(q, k, omega)