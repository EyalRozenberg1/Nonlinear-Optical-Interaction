clc
clear all
close all
%% Chosse which tolerance calculation 
AcceptTemperature   = 1;
AcceptBandwidth     = 1;

addpath ..\LBO

Lambda.in   = 1064.2e-9;    % Wavelength [m]
Lambda.out  = 1064.2e-9/2;  % Wavelength [m]
Pol.in      = 'z';
Pol.out     = 'y';
Temperature = 148.9009;
L           = 5e-2;         % Crystal Length

%% Temperature tolerance
syms T
% n2- n1
DeltaN(T) = n_lbo_Kato(Pol.out, Lambda.out, T) - n_lbo_Kato(Pol.in, Lambda.in, T);
Dn_DT     = diff(DeltaN,T);
Dn_DT_res = subs(Dn_DT, T , Temperature);

% Quasi PM SHG Tuning and Tolerance eq. (44)

DeltaT = 0.4429*Lambda.in/L * abs(Dn_DT_res).^-1;
fprintf('temperature tolerance is:  %f [Deg]\n',double(DeltaT));

%% Wavelength tolerance
syms lam
n.in(lam)   = n_lbo_Kato(Pol.in , lam, Temperature);
n.out(lam)  = n_lbo_Kato(Pol.out, lam, Temperature);

% n2 - n1
n2_minus_n1 = n_lbo_Kato(Pol.out, Lambda.out, Temperature) - n_lbo_Kato(Pol.in , Lambda.in, Temperature);

Dn1_Dlam = diff(n.in,lam);
Dn2_Dlam = diff(n.out,lam);

Dn1_Dlam_res = subs(Dn1_Dlam,lam,Lambda.in);
Dn2_Dlam_res = subs(Dn2_Dlam,lam,Lambda.out);

% Quasi PM SHG Tuning and Tolerance eq. (30)
DeltaLambda  = 0.4429*Lambda.in/L * abs(n2_minus_n1/Lambda.in + Dn1_Dlam_res - 0.5*Dn2_Dlam_res).^-1;
fprintf('BandWidth Tolerance is:    %f [um]\n',double(DeltaLambda./1e-6));
