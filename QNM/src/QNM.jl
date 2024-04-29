module QNM
export schroedinger!, bc!, xspan, uBoundaryNumerical, uHorizonNumerical, shear_mode_eq!, boundary_condition!, phiHorizonExpansion, dphiHorizonExpansion

#---------------------------------------------------------#
# Schrödinger equation #

# Wood-Saxon Potential
const V0 = 50 # MeV
const A = 40 # Massenzahl von Calcium
const R = 1.25 * A^(1 / 3) # fm
const a = 0.5 # fm

function V(x)
	return -V0 / (1 + exp((abs(x) - R) / a))
end

const ħ = 197 # Planck constant in MeV fm/c
const m = 939 # Nucleon mass in MeV/c^2
const xspan = (-15, 15) # fm

function schroedinger!(du, u, p, x)
	ψ = u[1]
	dψ = u[2]
	E = u[3]
	du[1] = dψ
	du[2] = 2 * m / ħ^2 * (V(x) - E) * ψ
	du[3] = 0
end

function bc!(residual, sol, p, x)
	residual[1] = sol(xspan[1])[1] # ψ(x=-15) = 0
	residual[2] = sol(xspan[end])[1] # ψ(x=15) = 0
end

#---------------------------------------------------------#
# Quasinormal modes #
# TODO test big numbers
epsilon = 1.0/1000000000.0
uBoundaryNumerical = epsilon
uHorizonNumerical = 1.0-1.0/10.0

f(u, qt) = 1 - (1+qt^2) * u^2 + qt^2 * u^3 # blackening_factor
df(u, qt) = -2 * (1+qt^2) * u + 3 * qt^2 * u^2 # derivative of blackening_factor

function shear_mode_eq!(du, u, p, t)
    # u is current state variable, du is the derivative of u at time t, t is current time
    c0, qt, kk = p # p is a vector of parameters

    phiReal, phiImag, dphiReal, dphiImag, omegaReal, omegaImag = u

    phi = phiReal + im * phiImag
    dphi = dphiReal + im * dphiImag
    omega = omegaReal + im * omegaImag

    du[1] = dphiReal
    du[2] = dphiImag
    eq = (
            -((omega^2 - kk^2*(-1 + t)*(-1 + t*(-1 + qt^2*t)))* phi)/(4*(-1 + t)^2*t*(1 + t - qt^2*t^2)^2) 
            -((-1 + t*(-1 + qt^2*t))*(-1 + t^2*(-1 + qt^2*(-1 + 2*t)))* dphi)/((-1 + t)*t*(1 + t - qt^2*t^2)^2)
        ) * (1 + t - qt^2*t^2)^2/((-1 + t*(-1 + qt^2*t))^2)
    du[3] = real(eq)
    du[4] = imag(eq)

    # Boundary condtions for omega
    du[5] = 0
    du[6] = 0
end

# boundary conditions
function boundary_condition!(residual, u, p, t)
    c0, qt, kk = p
    vars = [c0, qt, u[end][5] + im * u[end][6], kk]
    residual[1] = real(phiHorizonExpansion(uHorizonNumerical, vars...)) - u[end][1] # solution at horizon should be the horizon expansion
    residual[2] = imag(phiHorizonExpansion(uHorizonNumerical, vars...)) - u[end][2] # solution at horizon should be the horizon expansion
    residual[3] = real(dphiHorizonExpansion(uHorizonNumerical, vars...)) - u[end][3] # derivative at horizon should be the derivative of the horizon expansion
    residual[4] = imag(dphiHorizonExpansion(uHorizonNumerical, vars...)) - u[end][4] # derivative at horizon should be the derivative of the horizon expansion
end

#---------------------------------------------------------#
# Horizon expansion #

function phiHorizonExpansion(u, c0, qt, omega, kk)
	Sqrt = x -> sqrt(complex(x))
	I = im
	w = omega
	return (1 - u)^((im*w)/(2*(-2 + qt^2)))*
	 (c0 - (c0*(1 - u)*(kk^2*(-2 + qt^2)^2 + 
		  w*(2*Sqrt(-(-2 + qt^2)^2) - 4*w + 
					 qt^2*(2*Sqrt(-(-2 + qt^2)^2) + 5*w))))/(4*(-2 + 
		  qt^2)*(4 - 4*qt^2 + qt^4 - Sqrt(-(-2 + qt^2)^2)*w)) + 
		(c0*(1 - u)^2*(kk^4*(-2 + qt^2)^4 - 
		  2*kk^2*(-2 + qt^2)^2*(-36*qt^4 + 8*qt^6 + 
					 4*(-4 + Sqrt(-(-2 + qt^2)^2)*w + w^2) - 
			 qt^2*(-48 + 8*Sqrt(-(-2 + qt^2)^2)*w + 5*w^2)) + 
		  w*(-4*qt^8*(4*Sqrt(-(-2 + qt^2)^2) + 23*w) + 
			 4*qt^6*(4*Sqrt(-(-2 + qt^2)^2) + 103*w) + 
			 4*(16*Sqrt(-(-2 + qt^2)^2) - 8*w + 
				7*Sqrt(-(-2 + qt^2)^2)*w^2 + 4*w^3) - 
			 4*qt^2*(64*Sqrt(-(-2 + qt^2)^2) - 52*w + 
				21*Sqrt(-(-2 + qt^2)^2)*w^2 + 10*w^3) + 
					 
			 qt^4*(144*Sqrt(-(-2 + qt^2)^2) - 552*w + 
				80*Sqrt(-(-2 + qt^2)^2)*w^2 + 25*w^3))))/
		  (32*(-2 + qt^2)^2*(4 - 4*qt^2 + qt^4 - 
		  Sqrt(-(-2 + qt^2)^2)*w)*(8 - 8*qt^2 + 2*qt^4 - 
				Sqrt(-(-2 + qt^2)^2)*w)))
  end
  
  function dphiHorizonExpansion(u, c0, qt, omega, kk)
	  Sqrt = x -> sqrt(complex(x))
	  I = im
	  w = omega
	  return (1 - u)^((im*w)/(2*(-2 + qt^2)))*
	  ((c0*(kk^2*(-2 + qt^2)^2 + 
		  w*(2*Sqrt(-(-2 + qt^2)^2) - 4*w + 
			 qt^2*(2*Sqrt(-(-2 + qt^2)^2) + 5*w))))/
		   (4*(-2 + qt^2)*(4 - 4*qt^2 + qt^4 - 
		  Sqrt(-(-2 + qt^2)^2)*w)) - 
		 (c0*(1 - u)*(kk^4*(-2 + qt^2)^4 - 
		  2*kk^2*(-2 + qt^2)^2*(-36*qt^4 + 8*qt^6 + 
					  4*(-4 + Sqrt(-(-2 + qt^2)^2)*w + w^2) - 
			 qt^2*(-48 + 8*Sqrt(-(-2 + qt^2)^2)*w + 5*w^2)) + 
				 
		  w*(-4*qt^8*(4*Sqrt(-(-2 + qt^2)^2) + 23*w) + 
			 4*qt^6*(4*Sqrt(-(-2 + qt^2)^2) + 103*w) + 
					  
			 4*(16*Sqrt(-(-2 + qt^2)^2) - 8*w + 
				7*Sqrt(-(-2 + qt^2)^2)*w^2 + 4*w^3) - 
					  
			 4*qt^2*(64*Sqrt(-(-2 + qt^2)^2) - 52*w + 
				21*Sqrt(-(-2 + qt^2)^2)*w^2 + 10*w^3) + 
					  
			 qt^4*(144*Sqrt(-(-2 + qt^2)^2) - 552*w + 
				80*Sqrt(-(-2 + qt^2)^2)*w^2 + 25*w^3))))/
		   (16*(-2 + qt^2)^2*(4 - 4*qt^2 + qt^4 - 
		  Sqrt(-(-2 + qt^2)^2)*w)*(8 - 8*qt^2 + 2*qt^4 - 
				 Sqrt(-(-2 + qt^2)^2)*w))) - (1/(2*(-2 + 
		 qt^2)))*(I*(1 - u)^(-1 + (I*w)/(2*(-2 + qt^2)))*w*
		 (c0 - (c0*(1 - u)*(kk^2*(-2 + qt^2)^2 + 
			w*(2*Sqrt(-(-2 + qt^2)^2) - 4*w + 
						 
			   qt^2*(2*Sqrt(-(-2 + qt^2)^2) + 5*w))))/(4*(-2 + 
			qt^2)*(4 - 4*qt^2 + qt^4 - Sqrt(-(-2 + qt^2)^2)*w)) + 
			(c0*(1 - u)^2*(kk^4*(-2 + qt^2)^4 - 
			2*kk^2*(-2 + qt^2)^2*(-36*qt^4 + 8*qt^6 + 
						 4*(-4 + Sqrt(-(-2 + qt^2)^2)*w + w^2) - 
			   qt^2*(-48 + 8*Sqrt(-(-2 + qt^2)^2)*w + 5*w^2)) + 
					
			w*(-4*qt^8*(4*Sqrt(-(-2 + qt^2)^2) + 23*w) + 
			   4*qt^6*(4*Sqrt(-(-2 + qt^2)^2) + 103*w) + 
						 
			   4*(16*Sqrt(-(-2 + qt^2)^2) - 8*w + 
				  7*Sqrt(-(-2 + qt^2)^2)*w^2 + 4*w^3) - 
						 
			   4*qt^2*(64*Sqrt(-(-2 + qt^2)^2) - 52*w + 
				  21*Sqrt(-(-2 + qt^2)^2)*w^2 + 10*w^3) + 
			   qt^4*(144*Sqrt(-(-2 + qt^2)^2) - 552*w + 
				  80*Sqrt(-(-2 + qt^2)^2)*w^2 + 25*w^3))))/
			  (32*(-2 + qt^2)^2*(4 - 4*qt^2 + qt^4 - 
			Sqrt(-(-2 + qt^2)^2)*w)*(8 - 8*qt^2 + 2*qt^4 - Sqrt(-(-2 + qt^2)^2)*w))))
  end
  

end