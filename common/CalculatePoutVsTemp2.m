function [figNum, MaxTempValue] = CalculatePoutVsTemp2(PinPeak, Undepleted, T, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, TempSamples, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples)
% Calculating the BW of the conversion


A_from_I    = @(intensity,n)       sqrt(intensity/(n*eps0*c));
P_from_A    = @(A_r,r,RefIdx)      2*trapz(r, 2*pi*r.*(2*RefIdx*eps0*c*abs(A_r).^2)); % For a gaussian beam A is A(r)

TempVec = T.pm-2:0.1:T.pm+2;

% h           = waitbar(0,'calculating P2\omega vs T');
P2w= zeros(1,TempSamples);  % Memmory Allocation

for temp=1:TempSamples
    tempGrad = zeros(1,NumOfPoints) + TempVec(temp);

%     waitbar(temp / TempSamples);

    % Create Delta K
    [DeltaK, K, n, Omega] = DeltaK_Creator(tempGrad, InteractionType, Lambda, Pol, refIdx, k, w);

    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF2(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples,0);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AoutPower2  = A.out.*conj(A.out); AoutPower2  = AoutPower2(:,1).';
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        P2w(temp) = IOut(end);
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            disp('ERROR: numerical issue while calculating Pout Vs Temperature. Exit code')
            return;
        else
            P2w(temp) = P.out(end);
        end
    end
    fprintf('Finished calculation for sample number %d out of %d\n',temp,TempSamples);
end

% close(h);
% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, 0, 0, Lambda, TempVec, 0, P2w, 0, I, 0, InteractionType, Print, NumOfPoints, 0, 0, 0, 0);

% find the temperature for maximum efficiency
if(PlaneGauss_)
    P_G = P2w/10^4; % convertion from W/m^2 to W/cm^2
else
    P_G = P2w;
end
[~,S] = max(P_G);
MaxTempValue = TempVec(S);
end

