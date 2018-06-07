function [figNum, Tmax, Tmin, Tm] = CalculatePoutVsTm(figNum, GradType, T, TmSamples, CrystalPropAxis, NumOfPoints, L, InteractionType, Lambda, Pol, refIdx, k, w, PinPeak, Undepleted, w0, PlaneGauss_, I, dx_prop, c, eps0, deff, A_from_I, Kappa, Print, P_from_A, samples, xi)

% Tstep = linspace(-10,10,TmSamples);
DeltaBW     = 6;
Tstep       = -0.5*DeltaBW:DeltaBW/TmSamples:0.5*DeltaBW-DeltaBW/TmSamples; % Temperature

Efficiency          = zeros(1,TmSamples);  % Memmory Allocation

for Tm=1:TmSamples
    t.pm  = T.pm  + Tstep(Tm);
    t.max = T.max + Tstep(Tm);
    t.min = T.min + Tstep(Tm);
    
    [TempGrad]      = TemperatureGradient(t, L, CrystalPropAxis, GradType, NumOfPoints, xi);

    % Create Delta K
    [DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad, InteractionType, Lambda, Pol, refIdx, k, w);
    
    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples,0);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out);
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        Efficiency(Tm) = IOut(end)/I.in1;
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            disp('ERROR: numerical issue while calculating Pout Vs Temperature. Exit code')
            return;
        else
            Efficiency(Tm) = P.out(end)/P.in1(1);
        end
    end
    fprintf('Finished calculation for sample number %d out of %d\n',Tm,TmSamples);
end

% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, 0, Efficiency, Lambda, T.pm + Tstep, 0, 0, 0, I, 0, InteractionType, Print, NumOfPoints, 0, 0, 0, 0);

% find the temperature for maximum efficiency
[~,S] = max(Efficiency);
Tmax = T.max + Tstep(S);
Tmin = T.min + Tstep(S);
Tm   = T.pm + Tstep(S);
end

