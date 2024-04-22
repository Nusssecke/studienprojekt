# This copied straight from the mathematica notebook.
# Order 14
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