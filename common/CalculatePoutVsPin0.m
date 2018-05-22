function [figNum] = CalculatePoutVsPin0(PinPeak, Undepleted,w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, NumSamples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples)
% Calculating the BW of the conversion

%h       = waitbar(0,'calculating P2\omega vs p\omega(start)');

IinVec  = linspace(0,I.in1, NumSamples);
P2w     = zeros(1,NumSamples);  % Memmory Allocation
Pw      = zeros(1,NumSamples);  % Memmory Allocation

% Create Delta K
[DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad, InteractionType, Lambda, Pol, refIdx, k, w);
    %DeltaK = zeros(size(DeltaK)); % Perfect PM
for p=1:NumSamples
tic    
    i.in1 = IinVec(p);
    i.out = 0;
%     waitbar(p / NumSamples);

    % plane wave propagation using Split Step Fourier 
    [A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, i, InteractionType, deff, A_from_I, Kappa, P_from_A, samples,0);

    % intensity at the end of the crystal %
    if(PlaneGauss_)
        AinPower2   = A.in1.*conj(A.in1); AinPower2  = AinPower2(:,1).';
        Iin1        = 2*eps0*c*n.in1.*AinPower2;
        AoutPower2  = A.out.*conj(A.out); AoutPower2  = AoutPower2(:,1).';
        IOut        = 2*eps0*c*n.out.*AoutPower2;
        P2w(p)      = IOut(end);
        Pw(p)       = Iin1(1);
    else
        if(isnan(P.out(end)) || isnan(P.in1(end)))
            disp('ERROR: numerical issue while calculating Pout Vs Pin(start). Exit code')
            return;
        else
            P2w(p)  = P.out(end);
            Pw(p)   = P.in1(1);
        end
    end
    fprintf('Finished calculation for sample number %d out of %d\n', p ,NumSamples);
toc    
end

if(PlaneGauss_)
    Presult.in1 = Pw;
    Presult.out = P2w;
else
    Presult.in1 = Pw;
    Presult.out = P2w;   
end

% close(h);
% Print results
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, 0, 0, 0, Lambda, 0, 0, Presult, 0, I, 0, InteractionType, Print, NumOfPoints, 0, 0, 0, 0);
end

