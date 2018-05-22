function [T_mat, CEM_mat, dt] = solve_heat_equation_implicit_ADI(x, y, z, T, Qt, t_h, t_c, td, w, Cb, Ct, rho, rho_b, T_amb)

% This function solves three-dimensional heat equation in a homogeneous 
% medium using Alternating Direction Implicit (ADI) method.
%
% Inputs:
% 
% x = grid vector in x-direction (m)
% y = grid vector in y-direction (m)
% z = grid vector in z-direction (m)
% T = initial temperature matrix (C)
% Qt = heating rate matrix (K/s) (set to 0 if none)
% t_h = heating duration (s)
% t_c = cooling duration (s)
% td = thermal diffusivity (m^2/s)
% w = blood perfusion rate (volumetric flow rate per unit volume) (1/s) (set to 0 if none)
% Cb = specific heat capacity of blood (J/kg/K) (set to 0 if no perfusion)
% Ct = specific heat capacity of tissue (J/kg/K)
% rho = mass density of tissue (kg/m^3)
% rho_b = mass density of blood (kg/m^3) (set to 0 if no perfusion)
% T_amb = ambient temperature for blood perfusion (C) (set to 0 if no perfusion)
%
% Outputs:
%
% T_mat = Temperature matrix (C)
% CEM_mat = Thermal dose matrix (Cumulative Equivalent Minutes at 43 C)
% dt = time step (s)
%
% Visa Suomi
% Turku University Hospital
% 
% November 2017

%% compute parameters and initialise variables

% determine time step (s)

dt = min([0.1, (t_h+t_c)/10]);

% determine number of integration steps at each stage

N_h = round(t_h/dt);
N_c = round(t_c/dt);
N = N_h+N_c;                        % total number of time steps

% create heating on/off -vector

if(N_h~=0)
    initial = ones(1,N_h);
else
    initial = [];
end

if(N_c~=0)
    cooloff = zeros(1,N_c);
else
    cooloff = [];
end

input = [initial cooloff];

% determine grid resolution (m)

dx = abs(abs(x(1))-abs(x(2)));
dy = abs(abs(y(1))-abs(y(2)));
dz = abs(abs(z(1))-abs(z(2)));

% set constants

r1 = td*dt/(dx^2);
r2 = td*dt/(dy^2);
r3 = td*dt/(dz^2);

% set LHS matrix coefficients

a1 = -r1/2;
b1 = 1 + r1;
c1 = -r1/2;

a2 = -r2/2;
b2 = 1 + r2;
c2 = -r2/2;

a3 = -r3/2;
b3 = 1 + r3;
c3 = -r3/2;

%% Solve heat equation using ADI

T_mat = NaN([size(T) N]);       % matrix for storing temperature (C)
CEM_mat = NaN([size(T) N]);     % matrix for storing thermal dose (CEM43)

% start iteration over time

fprintf('Calculation started\n')
h = waitbar(0, 'Solving heat equation...');
tic
for n = 1:N
    
    % Execute x-direction
    
    T_s = T;                                % initialise temporary matrix
    
    for j = 2:length(y)-1
    for k = 2:length(z)-1
        
    % calculate RHS matrix
        
    d = NaN(length(x)-2,1);
        
    for i = 2:length(x)-1
        d(i-1) = r1/2 * T(i-1,j,k) + r1/2 * T(i+1,j,k) ...
                + r2 * T(i,j-1,k) + r2 * T(i,j+1,k) ...
                + r3 * T(i,j,k-1) + r3 * T(i,j,k+1) ...
                + (1 - r1 - 2*r2 - 2*r3) * T(i,j,k);
    end
    
    d(1) = d(1) - a1 * T(1,j,k);
    d(end) = d(end) - c1 * T(end,j,k);
        
    % Solve using tridiagonal matrix algorithm (Thomas algorithm)
    
    u = thomas_algorithm(a1, b1, c1, d);

    % fill solved u-vector into new arrays
    
    T_s(2:length(x)-1, j, k) = u;
    
    clear d u                               % clear temporary variables
    
    end
    end
    
    clear i j k                             % clear temporary variables
    
    % Execute y-direction
    
    T_ss = T;                               % initialise temporary matrix
    
    for i = 2:length(x)-1
    for k = 2:length(z)-1
        
    % calculate RHS matrix
        
    d = NaN(length(y)-2,1);
        
    for j = 2:length(y)-1
        d(j-1) = r1/2 * T(i-1,j,k) + r1/2 * T(i+1,j,k) ...
                + r1/2 * T_s(i-1,j,k) + r1/2 * T_s(i+1,j,k) ...
                + r2/2 * T(i,j-1,k) + r2/2 * T(i,j+1,k) ...
                + r3 * T(i,j,k-1) + r3 * T(i,j,k+1) ...
                - r1 * T_s(i,j,k) ...
                + (1 - r1 - r2 - 2*r3) * T(i,j,k);
    end
    
    d(1) = d(1) - a2 * T(i,1,k);
    d(end) = d(end) - c2 * T(i,end,k);
        
    % Solve using tridiagonal matrix algorithm (Thomas algorithm)
    
    u = thomas_algorithm(a2, b2, c2, d);
            
    % fill solved u-vector into new arrays
    
    T_ss(i, 2:length(y)-1, k) = u;
    
    clear d u                               % clear temporary variables
    
    end
    end
    
    clear i j k                             % clear temporary variables
    
    % Execute z-direction
    
    T_new = T;                              % initialise final matrix
    
    for i = 2:length(x)-1
    for j = 2:length(y)-1
        
    % calculate RHS matrix
        
    d = NaN(length(z)-2,1);
        
    for k = 2:length(z)-1
        d(k-1) = r1/2 * T(i-1,j,k) + r1/2 * T(i+1,j,k) ...
                + r1/2 * T_s(i-1,j,k) + r1/2 * T_s(i+1,j,k) ...
                + r2/2 * T(i,j-1,k) + r2/2 * T(i,j+1,k) ...
                + r2/2 * T_ss(i,j-1,k) + r2/2 * T_ss(i,j+1,k) ...
                + r3/2 * T(i,j,k-1) + r3/2 * T(i,j,k+1) ...
                - r1 * T_s(i,j,k) - r2 * T_ss(i,j,k) ...
                + (1 - r1 - r2 - r3) * T(i,j,k);
    end
    
    d(1) = d(1) - a3 * T(i,j,1);
    d(end) = d(end) - c3 * T(i,j,end);
        
    % solve using tridiagonal matrix algorithm (Thomas algorithm)
    
    u = thomas_algorithm(a3, b3, c3, d);
        
    % fill solved u-vector into new arrays
    
    T_new(i, j, 2:length(z)-1) = u;
    
    clear d u                               % clear temporary variables
    
    end
    end
    
    % update temperature
    
    Pt = rho_b*w*Cb/(rho*Ct)*(T_amb-T_new);
    T_new = T_new + dt*(Pt + input(n)*Qt);
    
    % calculate thermal dose
    
    Tind1 = double(T_new < 43);
    Tind2 = double(T_new >= 43);
    
    if max(Tind1(:)) == 1
        if n == 1
            CEM_mat(:,:,:,n) = dt*0.25.^(43-Tind1.*T_new);
        else
            CEM_mat(:,:,:,n) = CEM_mat(:,:,:,n-1) + dt*0.25.^(43-Tind1.*T_new);
        end
    end
    if max(Tind2(:)) == 1
        if n == 1
            CEM_mat(:,:,:,n) = dt*0.50.^(43-Tind2.*T_new);
        else
            CEM_mat(:,:,:,n) = CEM_mat(:,:,:,n-1) + dt*0.50.^(43-Tind2.*T_new);
        end
    end
    
    T_mat(:,:,:,n) = T_new;                 % update temperature matrix
    T = T_new;                              % initialise T for next iteration
    
    clear i j k T_s T_ss T_new Tind1 Tind2  % clear temporary variables
    
    % print progress  

    h = waitbar(n / N);
    
end
delete(h)
fprintf('Calculation finished\n')
toc