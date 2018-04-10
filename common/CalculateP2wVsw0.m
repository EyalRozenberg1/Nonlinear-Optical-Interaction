function [figNum, MaxWaistValue] = CalculateP2wVsw0(PinPeak, Undepleted, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, w0Samples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, I_from_A, Kappa, Print, P_from_A, samples)

% Gauss adiabatic depleted
BeamWidthVec_in1  = linspace(0.5*w0.in1,1.5*w0.in1,w0Samples); % Wavelength [m]
% BeamWidthVec_in1  = linspace(80e-6,160e-6,w0Samples); % Wavelength [m]

% h           = waitbar(0,'calculating P2w Vs w0');
Efficiency  = zeros(1,w0Samples);  % Memmory Allocation

% Create Delta K
[ DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad,InteractionType, Lambda, Pol, refIdx, k, w);
for beamwidth=1:w0Samples
tic
    w_0.in1        = BeamWidthVec_in1(beamwidth);
    w_0.out        = w_0.in1/sqrt(2);
    I.in1          = PinPeak/(pi*w_0.in1.^2);
%     waitbar(beamwidth / w0Samples);
    
    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, Lambda, w_0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out); AoutPower2  = AoutPower2(:,1).';
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        Efficiency(beamwidth) = IOut(end)/I.in1;
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            error('ERROR: numerical issue while calculating P2w Vs w0. Exit code');
        else
            Efficiency(beamwidth) = P.out(end)/P.in1(1);
        end
    end
    fprintf('Finished "P2w Vs w0" calculation for sample number %d out of %d\n',beamwidth,w0Samples);
toc
end
% close(h);

% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, 0, Efficiency, Lambda, 0, 0, 0,   0, I, 0, InteractionType, Print, NumOfPoints, 0, I_from_A, BeamWidthVec_in1, 0);

% find the waist for maximum efficiency

[~,S] = max(Efficiency);
MaxWaistValue.in1 = BeamWidthVec_in1(S);
MaxWaistValue.out = MaxWaistValue.in1/sqrt(2);

end

