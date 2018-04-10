function [ DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad,InteractionType, Lambda, Pol, refIdx, k, w)
% Creating:
% DeltaK & Kin1, Kin2, Kout
% In accordance to the type of frequency generation

switch InteractionType
    case 'Type1 SHG'
    n.in1 = refIdx(Pol.in1, Lambda.in1, TempGrad);
    n.in2 = n.in1;
    n.out = refIdx(Pol.out, Lambda.out, TempGrad);

    K.in1 = k(n.in1, Lambda.in1);
    K.in2 = K.in1;
    K.out = k(n.out, Lambda.out);

    Omega.in1 = w(Lambda.in1);
    Omega.in2 = Omega.in1;
    Omega.out = w(Lambda.out);

    % DeltaK = 2kin - kout
    DeltaK = 2*K.in1 - K.out;
%     DeltaK = DeltaK_CreatorGhosh(TempGrad, Lambda.in1);
    
end
            
end

