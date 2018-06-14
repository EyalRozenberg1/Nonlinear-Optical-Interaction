function [T, I, deff, Pol, Process, Lambda,...
    Walkoffangle, Orientation, L, dx_prop,...
    NumOfPoints, CrystalPropAxis, c, eps0,...
    k, w, refIdx, A_from_I, I_from_A,...
    Kappa, w0, P_from_A, samples, PinPeak,...
    Pulse_wdt, f, xi] = SetVal( InteractionType, Type )

% Setting all nedded values for the simulated interaction

% Crystal Input Parameters
L               = 5e-2;              % Crystal length [m] Default: 5e-2
dx_prop         = 5e-6;              % Split Step Fourier   Default: 5e-7
NumOfPoints     = round(L/dx_prop);  % round(L/dx_prop-1)
CrystalPropAxis = linspace(0,L,NumOfPoints);

% Gaussian Wave parameters
xi               = 5.1315; % {1, 2.84, 3.317, 5.1315}

% Intensity beam waist [m] {Alpha-Las:120e-6, Luce:49e-6 }
w0.I        = 72.639e-6/sqrt(2)/sqrt(xi);%49e-6/sqrt(2); % B&K for xi=2.84 -> w0.I=43.170e-06/sqrt(2)
w0.in1      = sqrt(2)*w0.I;         % fundamental  [m] - under calculation considering Lambda in   B&K {1/sqrt(2.84)*72.75e-6 ~= 43.169e-06}
w0.out      = w0.in1/sqrt(2);       % 2nd harmonic [m] - under calculation considering Lambda out  B&K {1/sqrt(2.84)*51.44e-6 ~= 30.524e-05}
samples     = 200; % 400            % Hankel number of sumples definition
Attenuation = 1;                  % 0.56; 1/sqrt(2*pi); energy attenuation due to pusle gaussian shape


% General Parameters
c           = 299792458;               % m/s Speed of light
eps0        = 8.854187817e-12;         % F/m
k           = @(n, lambda)              n.*2*pi/lambda;
w           = @(lambda)                 2*pi*c/lambda;
refIdx      = @(pol, lambda, temp)      n_lbo_Kato(pol, lambda, temp);
% refIdx      = @(pol, lambda, temp)      n_lbo_Ghosh(pol, lambda, temp);
Kappa       = @(d_eff,w_,k_,r)          1i*w_^2*d_eff./(k_(r)*c^2)*Attenuation;

% according to boyd pg. 98 equation 2.7.39:
A_from_I    = @(intensity,n)       sqrt(intensity./(2*n*eps0*c));
I_from_A    = @(RefIdx,AvPower2)   2*eps0*c*RefIdx.*AvPower2;
% according to boyd pg. 118 (end of the page):
P_from_A    = @(A_r,r,n)           trapz(r, 2*pi*r.*(2*n*eps0*c*abs(A_r).^2)); % For a gaussian beam A is A(r)

switch InteractionType
    case 'Type1 SHG'
    if strcmp(Type,'non critical')
        Lambda.in1  = 1064.2e-9;   % Wavelength [m]
        Lambda.in2  = Lambda.in1;
        Lambda.out  = Lambda.in1/2;
        PinAvg      = 300e-3;  % [W]   Luce: 940e-3[W]   Alpha Las: 100e-3[W]   fiber laser: 6[W]
        f           = 10e3;    % [Hz]  Luce: 10K[Hz]     Alpha Las: 100[Hz]     fiber laser: 20K[Hz]
        Pulse_wdt   = 5e-9;    % [s]   Luce: 5n[sec]     Alpha Las: 1.1n[sec]   fiber laser: 10n[sec]
        PinPeak     = PinAvg/(f*Pulse_wdt);
        I.in1       = PinPeak/(pi*w0.I^2);% Input  Intensity 1[W/m^2]
        I.in2       = I.in1;        % Input  Intensity 2[W/m^2]
        I.out       = 0;            % Output Intensity  [W/m^2]
        Pol.in1     = 'z';
        Pol.in2     = 'z';
        Pol.out     = 'y';
        Process     = 'oo->e';
        deff        = 1.05e-12;     % d32[m/V] Raicol:0.85e-12 Coherent:1.2e-12 Laser Components:-0.98e-12 Extma:1.05e-12+-0.09e-12
        Walkoffangle= 0;            % [mrad]
        Orientation.theta = 90;     % [Deg]
        Orientation.phi   = 0;      % [Deg]
        T.pm              = 149.3709;% 148.900862: for Dk ~=0, 149.3409: for optimum efficiency
        T.min             = T.pm - 0.3492; %20 start of the crystal; xi=3.317, -0.24
        T.max             = T.pm + 0.3492; %20 end of the crystal; xi=3.317, +0.24
%         xi=5,     T.min=-0.28, T.max=+0.28, T.pm=149.265
%         xi=3.317, T.min=-0.24, T.max=+0.24, T.pm=149.255
%         xi=2.84,  T.min=-0.23, T.max=+0.23, T.pm=149.247
%         T.min             = T.pm - 0.4594;     % LBO Crystal low gradient Temperature.
%         T.max             = T.pm + 0.4594;     % LBO Crystal high gradient Temperature.
        
        % Adiabatic:
        % Plane - I.in1 = 8.8174e+12[W/m^2] (A.las)
        %         T.min = T.pm + 5.655;
        %         T.max = T.pm - 5.655;
        % Gauss - P.in1 Peak = 1M[W] (A.las)
        %         T.min = T.pm + 9.5489;
        %         T.max = T.pm - 9.5489;
        
        
        % Gauss - Good results for Const 149.3114 Pin=30KW
        % Gauss - Good results for Const 149.4741
        % Gauss - Good results for Linear:       T.min=144       & T.max=153
        % Gauss - Good results for Apodization1: T.min=149.315   & 149.303
        
        % Plane - Good results for Const 148.900862 (the old was 148.9009) to get DelkaK ~= 0
        % Plane - Good results for Apodization1: T.min=140   & T.max=156.5
        % Plane - Good results for Linear:       T.min=143.5 & T.max=154 or
        %                                        T.min=140   & T.max=156.5
    end
end

xi = L/(2*pi*1.605*w0.in1^2/Lambda.in1);

end

