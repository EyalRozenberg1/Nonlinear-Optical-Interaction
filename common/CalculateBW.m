function [figNum] = CalculateBW(PinPeak, Undepleted, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, BWSamples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, I_from_A, Kappa, Print, P_from_A, samples, T)
% Calculating the BW of the conversion

% 1) adiabatic
% LambdaVec.in1 = linspace(Lambda.in1-15e-9,Lambda.in1+20e-9,BWSamples); % Wavelength [m]
% 1) adiabatic
% LambdaVec.in1 = linspace(Lambda.in1-25e-9,Lambda.in1+25e-9,BWSamples); % Wavelength [m]
DeltaBW = 150e-9;
LambdaVec.in1 = Lambda.in1-0.5*DeltaBW:DeltaBW/BWSamples:Lambda.in1+0.5*DeltaBW-DeltaBW/BWSamples; % Wavelength [m]



% h           = waitbar(0,'calculating BW');
Efficiency  = zeros(1,BWSamples);  % Memmory Allocation

% in case of a constant temperature gradient
DeltakL  = 0;%zeros(1,BWSamples);  % Memmory Allocation;

% Create constant temperature gradient for DeltaK*L axis
% [tempGradConst] = TemperatureGradient( T, CrystalPropAxis(end), CrystalPropAxis, 'Const', NumOfPoints);
% Crystal length
% L = CrystalPropAxis(end);

for lam=1:BWSamples
tic
    lam
    lambda.in1 = LambdaVec.in1(lam);
    lambda.in2 = lambda.in1;
    lambda.out = lambda.in1/2;
%     waitbar(lam / BWSamples);

    % Create Delta K
    [DeltaK, K, n, Omega]    = DeltaK_Creator(TempGrad,InteractionType, lambda, Pol, refIdx, k, w);
    
    % const Delta K for axis value only
%     [ DeltakConst, ~, ~, ~] = DeltaK_Creator(tempGradConst,InteractionType, lambda, Pol, refIdx, k, w);
%     DeltakL(lam) = DeltakConst(round(NumOfPoints/2))*L;

    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, lambda, w0,  NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples);
    
    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out);
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        Efficiency(lam) = IOut(end)/I.in1;
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            error('ERROR: numerical issue while calculating BW. Exit code');
        else
            Efficiency(lam) = P.out(end)/P.in1(1);
        end
    end
    fprintf('Finished BW calculation for sample number %d out of %d\n',lam,BWSamples);
toc
end

% close(h);
% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, DeltakL, Efficiency, LambdaVec, 0, 0, 0, 0, I, 0, InteractionType, Print, NumOfPoints, Print.BW, I_from_A, 0, 0);
end

