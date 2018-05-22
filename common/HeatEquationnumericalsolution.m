% This script compares the heat equation solvers with analytical 
% solution for radial heat source from "Kevin J. Parker, Effects of heat 
% conduction and sample size on ultrasonic absorption measurement, 
% JASA, 1985".
%
% Visa Suomi
% Turku University Hospital
%
% November 2017

%% initialise

close all;
clearvars;

%% define parameters

att = 0.1151;                           % attenuation coefficient (Np/cm)
rho = 1.05;                             % density tissue (g/cm^3)
rho_b = 1.05;                           % density blood (g/m^3)
c = 3.7630;                             % specific heat tissue (J/g/K)
I0 = 318;                               % maximum intensity (W/cm^2)
beta = 9e-2;                            % gaussian variance (cm^2)
tc = 0.51;                              % thermal conductivity (W/m/K)
dx = 0.1;                               % x-resolution (cm)
dy = 0.1;                               % y-resolution (cm)
dz = 0.1;                               % z-resolution (cm)
x = -1:dx:1;                            % x-vector (cm)
y = -1:dy:1;                            % y-vector (cm)
z = -1:dz:1;                            % z-vector (cm)
dur_heat = 0.9;                         % heating duration (s)
dur_cool = 9;                           % cooling duration (s)

%% create cylindar heating rate matrix

I = NaN(length(x), length(y));

for i = 1:length(x)
for j = 1:length(y)
    r = sqrt(x(i)^2 + y(j)^2);
    I(i,j) = I0*exp(-r.^2/beta);
end
end

I = repmat(I, 1, 1, length(z));         % intensity (W/cm^2)
Qt = 2*att*I/(rho*c);                   % heating rate (K/s)

% define boundary conditions

Qt(1,:,:) = 0;
Qt(end,:,:) = 0;
Qt(:,1,:) = 0;
Qt(:,end,:) = 0;
Qt(:,:,1) = 0;
Qt(:,:,end) = 0;

%% solutions over time at origin (r=0)

w = 0;                                  % perfusion rate (1/s)
Cb = 3617;                              % specific heat capacity of blood (J/kg/K)
Ct = c*1e3;                             % specific heat capacity of tissue (J/kg/K)
td = tc/(rho*1e3*Ct);                   % thermal diffusivity (m^2/s)
k = td/1e-4;                            % thermal diffusivity (cm^2/s)
T_amb = 0;                              % ambient temperature (C)
T_ini = T_amb*ones(size(Qt));           % initial temperature matrix (C)

% ADI method

[T_adi, ~, dt] = solve_heat_equation_implicit_ADI(x/1e2, y/1e2, z/1e2,...
    T_ini, Qt, dur_heat, dur_cool, td, w, Cb, Ct, rho*1e3, rho_b*1e3, T_amb);
T_adi_orig = squeeze(T_adi(median(1:length(x)), median(1:length(y)), median(1:length(z)), :))';

% analytical

t_vec_heat = dt:dt:dur_heat;
T_ana_orig = 2*att*I0*beta/(rho*c*4*k)*log(1 + 4*k*t_vec_heat/beta);

%% plotting

figure;
plot([0 t_vec_heat], [T_amb T_adi_orig(1:length(t_vec_heat))], '-.', 'LineWidth', 2); hold on;
plot([0 t_vec_heat], [T_amb T_ana_orig]); grid on;
xlabel('Time (s)'); ylabel('Temperature (C)');
legend('ADI', 'Analytical', 'Location', 'NW');
title('Temperature with time at origin');

%% solutions over radius during heating 
% (small heating duration only, see error in paper for the analytical solution)

heat_vec = dur_heat/3:dur_heat/3:dur_heat;          % heating times (s)
T_adi_heat = NaN(length(x), length(heat_vec));      % for storing values
T_ana_heat = NaN(length(x), length(heat_vec));

for i = 1:length(heat_vec)
    
    % ADI method
    
    T_adi_heat(:,i) = T_adi(:, median(1:length(y)), median(1:length(z)), round(heat_vec(i)/dt));
    
    % analytical
    
    T_ana_heat(:,i) = 2*att*I0*heat_vec(i)/(rho*c*(1+4*k*0/beta))*exp(-x.^2/(4*k*0+beta));
    
end

%% plotting

legendlabels = cell(2, length(heat_vec));
plot_col = {'r', 'g', 'b', 'k', 'c', 'm'};

figure;
for i = 1:length(heat_vec)
    
    plot(x, T_adi_heat(:, i), [plot_col{i} '-.'], 'LineWidth', 2); hold on;
    plot(x, T_ana_heat(:, i), [plot_col{i} '-']); hold on;
    legendlabels{1, i} = ['ADI ' num2str(heat_vec(i)) ' s']; 
    legendlabels{2, i} = ['Analytical ' num2str(heat_vec(i)) ' s'];

end

grid on;
xlabel('Radius (cm)'); ylabel('Temperature (C)');
title('Temperature profile over radius (heating)')
legend(legendlabels(:));
    
%% Analytical and implicit solution over radius during cooling

cool_vec = dur_cool/3:dur_cool/3:dur_cool;          % cooling times (s)

T_adi_cool = NaN(length(x), length(cool_vec));      % for storing values
T_ana_cool = NaN(length(x), length(cool_vec));

for i = 1:length(cool_vec)

    T_adi_cool(:,i) = T_adi(:, median(1:length(y)), median(1:length(z)), round((heat_vec(end) + cool_vec(i))/dt));
    T_ana_cool(:,i) = 2*att*I0*dur_heat/(rho*c*(1+4*k*cool_vec(i)/beta))*exp(-x.^2/(4*k*cool_vec(i)+beta));
    
end

%% plotting

figure;
for i = 1:length(cool_vec)
    
    plot(x, T_adi_cool(:, i), [plot_col{i} '-.'], 'LineWidth', 2); hold on;
    plot(x, T_ana_cool(:, i), [plot_col{i} '-']); hold on;
    legendlabels{1, i} = ['ADI ' num2str(cool_vec(i)) ' s'];
    legendlabels{2, i} = ['Analytical ' num2str(cool_vec(i)) ' s'];

end

grid on;
xlabel('Radius (cm)'); ylabel('Temperature (C)');
title('Temperature profile over radius (cooling)')
legend(legendlabels(:));