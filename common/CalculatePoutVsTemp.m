function [figNum, MaxTempValue] = CalculatePoutVsTemp(PinPeak, Undepleted, T, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, TempSamples, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples)
% Calculating the BW of the conversion

% TempVec = linspace(0.973*T.pm,1.02758*T.pm,TempSamples);
% TempVec = linspace(T.pm-3,T.pm+3,TempSamples);
% TempVec = T.pm-2:0.1:T.pm+2;
% TempVec   = linspace(T.pm-10,T.pm+10,TempSamples);
DeltaBW   = 6;
TempVec   = T.pm-0.5*DeltaBW:DeltaBW/TempSamples:T.pm+0.5*DeltaBW-DeltaBW/TempSamples; % Temperature
% TempVec = 146.9:0.1:151.9;

DeltaKL = zeros(size(TempVec)); % Memory allocation
L       = CrystalPropAxis(end);
% h           = waitbar(0,'calculating P2\omega vs T');
P2w         = zeros(1,TempSamples);  % Memmory Allocation
Efficiency  = zeros(1,TempSamples);  % Memmory Allocation
for temp=1:TempSamples
    tempGrad = zeros(NumOfPoints,1) + TempVec(temp);

%     waitbar(temp / TempSamples);

    % Create Delta K
    [DeltaK, K, n, Omega] = DeltaK_Creator(tempGrad, InteractionType, Lambda, Pol, refIdx, k, w);
    DeltaKL(temp) = DeltaK(1)*L;
    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out);
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        P2w(temp)   = IOut(end);
        Efficiency(temp) = IOut(end)/I.in1;
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            disp('ERROR: numerical issue while calculating Pout Vs Temperature. Exit code')
            return;
        else
            P2w(temp)        = P.out(end);
            Efficiency(temp) = P.out(end)/P.in1(1);
        end
    end
    fprintf('Finished calculation for sample number %d out of %d\n',temp,TempSamples);
end

% close(h);
% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, DeltaKL, Efficiency, Lambda, TempVec, 0, P2w, 0, I, 0, InteractionType, Print, NumOfPoints, 0, 0, 0, 0);

% find the temperature for maximum efficiency
if(PlaneGauss_)
    P_G = P2w/10^4; % convertion from W/m^2 to W/cm^2
else
    P_G = P2w;
end
[~,S] = max(P_G);
MaxTempValue = TempVec(S);
end

