using QNM, OrdinaryDiffEq, BoundaryValueDiffEq, Plots, DifferentialEquations, BenchmarkTools

# Rodas4P() not good enough
# Rosenbrock23() Quite good 2 decimal placement
# Rodas5() # fast
# Rodas5P()
# ESDIRK547L2SA2() very bad

# Test different solvers with computeQNM
solvers = [
    Euler(),
    Midpoint(),
    Heun(),
    Ralston(),
    RK4(),
    BS3(),
    BS5(),
    OwrenZen3(),
    OwrenZen4(),
    OwrenZen5(),
    DP5(),
    Tsit5(),
    RK065(),
    TanYam7(),
    DP8(),
    TsitPap8(),
    Feagin10(),
    Feagin12(),
    Feagin14(),
    MSRK5(),
    MSRK6(),
    Stepanov5(),
    SIR54(),
    Alshina2(),
    Alshina3(),
    Alshina6(),
    Vern6(),
    Vern7(),
    Vern8(),
    Vern9(),

    # Conservation laws?

    # Low storage methods
    ORK256(), # still starts using too much storage
    SSPRK53_2N12(),
    SSPRK53_2N2(),
    CarpenterKennedy2N54(),
    NDBLSRK124(),

    # Parallelized Explicit Extrapolation Methods
    AitkenNeville(),
    ExtrapolationMidpointDeuflhard(),

    # Explicit Multistep Methods
    AB3(),
    AB4(),
    AB5(),
    ABM54(),

    VCAB3(),

    # Stiff equations
    ImplicitEuler(),
    Cash4(),
    
    RadauIIA5(),

    Rodas5P(),

    Rosenbrock23(),
]