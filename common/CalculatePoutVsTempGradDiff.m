function [figNum, Tmax, Tmin] = CalculatePoutVsTempGradDiff(figNum, GradType, T, GradSamples, CrystalPropAxis, NumOfPoints, L, InteractionType, Lambda, Pol, refIdx, k, w, PinPeak, Undepleted, w0, PlaneGauss_, I, dx_prop, c, eps0, deff, A_from_I, Kappa, Print, P_from_A, samples)
% Calculating the BW of the conversion

% For adiabatic:
% DeltaT = linspace(-0.5*(T.max-T.min),0.5*(T.max-T.min),GradSamples);

% For Gouy Phase Compensation:
DeltaT = linspace(-(T.max-T.min),(T.max-T.min),GradSamples);

% h           = waitbar(0,'calculating P2\omega vs T');
Efficiency          = zeros(1,GradSamples);  % Memmory Allocation
DeltaKAxis          = zeros(1,GradSamples);  % Memmory Allocation

t.pm  = T.pm;
for diff=1:GradSamples
    t.max = T.pm + DeltaT(diff);
    t.min = T.pm - DeltaT(diff);
    
    [TempGrad]      = TemperatureGradient(t, L, CrystalPropAxis, GradType, NumOfPoints);

%     waitbar(temp / TempSamples);

    % Create Delta K
    [DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad, InteractionType, Lambda, Pol, refIdx, k, w);

    DeltaKAxis(diff) = (DeltaK(end) - DeltaK(1))*L;
    
    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out); AoutPower2  = AoutPower2(:,1).';
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        Efficiency(diff) = IOut(end)/I.in1;
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            disp('ERROR: numerical issue while calculating Pout Vs Temperature. Exit code')
            return;
        else
            Efficiency(diff) = P.out(end)/P.in1(1);
        end
    end
    fprintf('Finished calculation for sample number %d out of %d\n',diff,GradSamples);
end

% close(h);
% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, DeltaKAxis, Efficiency, Lambda, DeltaT, 0, 0, 0, I, 0, InteractionType, Print, NumOfPoints, 0, 0, 0, 0);

% find the temperature for maximum efficiency
[~,S] = max(Efficiency);
Tmax = T.pm + DeltaT(S);
Tmin = T.pm - DeltaT(S);
end

