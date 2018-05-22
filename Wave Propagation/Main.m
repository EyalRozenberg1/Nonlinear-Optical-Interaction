
% The code simulates a Efficient, Broadband, Robust frequency converstion
% using a LBO Crystal.

% The phase match will be achived by a temperature gradient
% Solution will include a adiabatic three-wave mixing solution.


% close all
% clc
addpath ..\common
addpath ..\LBO
addpath ..\Print
addpath ..\ExpResults

% Define required interaction process:
InteractionType = 'Type1 SHG';
% options - 
% Type1 SHG (available)
% Type2 SHG (not available)
% Type1 THG (not available)
% Type2 THG (not available)
Type = 'non critical';
% Options -
% non critical (available) - only for SHG
% critical     (not available)


% set simulation values according to interaction
A = 0; 
P = 0;
[T, I, deff, Pol, Process, Lambda,...
    Walkoffangle, Orientation, L, dx_prop,...
    NumOfPoints, CrystalPropAxis, c, eps0,...
    k, w, refIdx, A_from_I, I_from_A,...
    Kappa, w0, P_from_A, samples, PinPeak,...
    Pulse_wdt, f, xi] = SetVal( InteractionType, Type );



% What are the required prints
[ Print,figNum ,PlaneGauss_, GradType, Undepleted] = ReqPrints();


% first define the gradient function type:
[TempGrad]      = TemperatureGradient(T, L, CrystalPropAxis, GradType, NumOfPoints, samples);

if(Print.P2w_vs_Temp) % Only for constant temperature gradient - the tempgradient is set in the function
TempSamples            = 200; % 300               % Number of samples
[figNum, MaxTempValue] = CalculatePoutVsTemp (PinPeak, Undepleted, T, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, TempSamples, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples);
% [figNum, MaxTempValue] = CalculatePoutVsTemp2(PinPeak, Undepleted, T, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, TempSamples, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples);
Print.P2w_vs_Temp=0;
T.pm            = MaxTempValue;
if(strcmp(GradType,'Const'))
    [TempGrad]      = TemperatureGradient(T, L, CrystalPropAxis, GradType, NumOfPoints);
end
end


if(Print.P2w_vs_GradDiff)
GradSamples          = 200;               % Number of samples
[figNum, Tmax, Tmin] = CalculatePoutVsTempGradDiff(figNum, GradType, T, GradSamples, CrystalPropAxis, NumOfPoints, L, InteractionType, Lambda, Pol, refIdx, k, w, PinPeak, Undepleted, w0, PlaneGauss_, I, dx_prop, c, eps0, deff, A_from_I, Kappa, Print, P_from_A, samples);
Print.P2w_vs_GradDiff=0;
% T.max            = Tmax;
% T.min            = Tmin;
% [TempGrad]       = TemperatureGradient(T, L, CrystalPropAxis, GradType, NumOfPoints);
end


if(Print.P2wVsTm)
TmSamples          = 200;               % Number of samples
[figNum, Tmax, Tmin] = CalculatePoutVsTm(figNum, GradType, T, TmSamples, CrystalPropAxis, NumOfPoints, L, InteractionType, Lambda, Pol, refIdx, k, w, PinPeak, Undepleted, w0, PlaneGauss_, I, dx_prop, c, eps0, deff, A_from_I, Kappa, Print, P_from_A, samples);
Print.P2wVsTm=0;
% T.max            = Tmax;
% T.min            = Tmin;
% [TempGrad]       = TemperatureGradient(T, L, CrystalPropAxis, GradType, NumOfPoints);
end

if(Print.P2wVsw0)
w0Samples           = 40;                  % Number of samples
[figNum, MaxWaistValue] = CalculateP2wVsw0(PinPeak, Undepleted, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, w0Samples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, I_from_A, Kappa, Print, P_from_A, samples);
Print.P2wVsw0=0;
w0.in1 = MaxWaistValue.in1;
w0.out = MaxWaistValue.out;
w0.I   = w0.in1/sqrt(2); % intensity waist
I.in1  = PinPeak/(pi*w0.I.^2);
end

if(Print.P2w_vs_Pw)
NumSamples           = 20;%100                 % Number of samples
[figNum] = CalculatePoutVsPin0(PinPeak, Undepleted, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, NumSamples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, Kappa, Print, P_from_A, samples);
Print.P2w_vs_Pw=0;
end

if(Print.BW)
BWSamples           = 200;                  % Number of samples
[figNum] = CalculateBW(PinPeak, Undepleted, w0, PlaneGauss_, figNum, NumOfPoints, I, CrystalPropAxis, dx_prop, BWSamples, TempGrad, Lambda, Pol, refIdx, k, w, c, eps0, InteractionType, deff, A_from_I, I_from_A, Kappa, Print, P_from_A, samples, T);
Print.BW=0;
end

% Calculating results for given wavelength
if(Print.DeltaK || Print.Temperature || Print.RefIndex || Print.Amplitude || Print.Intensity || Print.NormST || Print.GouyPhase)
% Create Delta K
[DeltaK, K, n, Omega] = DeltaK_Creator(TempGrad,InteractionType, Lambda, Pol, refIdx, k, w);
end

% plane wave propagation using Split Step Fourier
if(Print.Amplitude || Print.Intensity || Print.NormST)
[A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples, Print.NormST);
end

% Print results
if(Print.DeltaK || Print.Temperature || Print.RefIndex || Print.Amplitude || Print.Intensity || Print.GouyPhase)
[figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, CrystalPropAxis, DeltaK, 0, Lambda, T, A, P, n, I, TempGrad, InteractionType, Print, NumOfPoints, Print.BW, I_from_A, 0, xi);
end

% Run normalized sationary states
if(Print.NormST)
[figNum] = sim_and_stationary_states_norm_units_for_JOSAB(CrystalPropAxis, A, Omega, K, DeltaK, deff, c, figNum);
end

if(Print.Exp.Results)
date      = '19_09_17';
% f         = 10*10^3;  % [Hz]   Lus: 10K[Hz]  Alpha Las: 100[Hz]
% Pulse_wdt = 5*10^-9; % [s]     Lus: 5n[sec]  Alpha Las: 1.1n[sec]
[figNum] = PrintExpResults(Print, figNum, date, f, Pulse_wdt);
end
