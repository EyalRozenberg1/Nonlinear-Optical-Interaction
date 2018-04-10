
% Based on Ghosh, Journal of Applied Physics, 78 (11) 1995

function [ DeltaK ] = DeltaK_CreatorGhosh(TempGrad, Lambda)
    
    Delta_n = -1.36245e-3 + 4.64015e-6*TempGrad + 5.84559e-8*TempGrad.^2 - 1.89162e-10*TempGrad.^3;
    DeltaK  = 2*2*pi/Lambda*( Delta_n );
                
end

