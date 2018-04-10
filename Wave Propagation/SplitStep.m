% Split step fourier solution for Plane wave propagation
% Generated 1/10/2017 By Eyal Rozenberg

function [P] = SplitStep(Undepleted, Lambda, w0, NumOfPoints, dz_prop, DeltaK, k_in1, k_out, Omega, n_in1, n_out, I, deff, samples)

% inputs:
% -------
% Undepleted:       pump depletion - 1 undepleted, 0 depleted
% Lambda:           Fundamental wave length [m]
% w0:               beam waists [m]
% NumOfPoints:      Number of propagation points (crystal length = NumOfPoints*dz_prop)
% dz_prop:          Propagation step [m]
% DeltaK:           wave vector mismatch - 2k1-k2 - vector lentgh NumOfPoints
% k_in1             Fundamental wavenumber - vector lentgh NumOfPoints
% k_out             SH wavenumber - vector
% Omega             FW frequency
% n_in1             FW refractive index - vector lentgh NumOfPoints
% n_out             SH refractive index - vector lentgh vector lentgh NumOfPoints
% I                 FW intensity [W/m^2]
% deff              nonlinear component
% samples           number of hankle transform samples

% Outpus:
% ---------
% P.in1             FW Power [W] - vector
% P.out             SH Power [W] - vector

c           = 299792458;               % m/s Speed of light
eps0        = 8.854187817e-12;         % F/m

P_from_A    = @(A_r,r,n)           trapz(r, 2*pi*r.*(2*n*eps0*c*abs(A_r).^2)); % For a gaussian beam A is A(r)
A_from_I    = @(intensity,n)       sqrt(intensity/(2*n*eps0*c));
Kappa       = @(d_eff,w_,k_,r)          1i*w_^2*d_eff/(k_(r)*c^2);

tic 
    x_max       = NumOfPoints * dz_prop;
    % Hankel definitions
    % Note: we can be defined some parameters as const because change in very little
    Kin1        = k_in1(1);

    zr.in1      = pi*n_in1(1)*w0^2/Lambda;     % [m] Reighley range with consider of Lambda in
    w_max.in1   = w0*sqrt(1+(x_max/zr.in1)^2);

    rmax.in1      = 2*w_max.in1;
    rmax.max      = rmax.in1;

    mat_H2        = hankel_matrix2(0, rmax.max, samples);
    r_            = mat_H2.r;                  % Radial co-ordinate vector
    fr_sq_        = (mat_H2.v).^2;             % frequency co-ordinate vector

    % Hankel Transform & Inverse Hankel Transform for out wave
    ht2         = @(f) mat_H2.JV .* (mat_H2.T * (f./mat_H2.JR)) /(2*pi);
    iht2        = @(F) mat_H2.JR .* (mat_H2.T * (F./mat_H2.JV)) *(2*pi);

    % propagation H same as Split Step
    H              = @(PropT, k)  exp(PropT).*exp(-1i*k*dz_prop);
    PropTerm       = @(lambda,n_) 1i*dz_prop*2*pi*sqrt((n_/lambda)^2 - fr_sq_);

    % according to boyd pg. 118 equations 2.10.5b-c:
    b_             = @(k,w_0) k*w_0^2;
    b.in1          = b_(Kin1, w0);

    x_waist        = x_max/2; % [m] waists location in the crystal
    Points_waist   = round(x_waist/dz_prop); % % waists point number in the propadation points
    chi_           = @(B) 2*(0-x_waist)/B;
    chi.in1        = chi_(b.in1);

    % Result for the start of the crystal
    % memory allocation
    P.in1 = complex(zeros(1,NumOfPoints));
    P.out = complex(zeros(1,NumOfPoints));

    % according to boyd pg. 117 equation 2.10.5a:
    Ain1   = A_from_I(I,n_in1(1))/(1+1i*chi.in1)*exp(-(r_.^2)/(w0^2*(1+1i*chi.in1)));
    Aout   = complex(zeros(samples,1));

    P.in1(1) = P_from_A(Ain1, r_, n_in1(1));
    P.out(1) = 0;

%                 for i=2:NumOfPoints changed
    for i=1:NumOfPoints

        % Prop term
        propterm.in1    = PropTerm(Lambda ,n_in1(i));
        propterm.out    = PropTerm(Lambda/2 ,n_out(i));

        propterm.in1( real(propterm.in1)>0 ) = -propterm.in1( real(propterm.in1)>0 );
        propterm.out( real(propterm.out)>0 ) = -propterm.out( real(propterm.out)>0 );

        Ew_in1 = ht2(Ain1);
        Ew_out = ht2(Aout);

        Hw.in1 = H(propterm.in1, k_in1(i));
        Hw.out = H(propterm.out, k_out(i));                    

        Ew_in1 = Ew_in1 .* Hw.in1;
        Ew_out = Ew_out .* Hw.out;

        Ain1_temp = iht2(Ew_in1);
        Aout_temp = iht2(Ew_out);

        % coupled equations, boyd pg. 98 equations 2.7.10-12:
        % General Definitions
        KappaIn  = Kappa(deff,Omega,k_in1,i);
        Kappaout = Kappa(deff,2*Omega,k_out,i);

%                     xi    = i * dz_prop; changed
        xi    = ( i - Points_waist ) * dz_prop;
        dAin1 = 2*KappaIn*Aout.*conj(Ain1)*exp( -1i*DeltaK(i)*xi );
        dAout = Kappaout*Ain1.^2*exp( 1i*DeltaK(i)*xi );

        if(Undepleted)
            Ain1 = Ain1;
        else
            Ain1 = Ain1_temp + dAin1*dz_prop;
        end
        Aout = Aout_temp + dAout*dz_prop;

        if(Undepleted)
            P.in1(i) = P.in1(1); % Undepleted Pump
        else
            P.in1(i) = P_from_A(Ain1, r_, n_in1(i));                        
        end
        P.out(i) = P_from_A(Aout, r_, n_out(i));

        if(isnan(P.out(i)) || isnan(P.in1(i)))
            error('ERROR: numerical issue while calculating BW. Exit code');
        end
        if(mod(i,round(NumOfPoints/10))==0 || i==NumOfPoints || i==2)
            disp(['Calc Split step. Done: ', num2str(round(100*i/NumOfPoints)),'%',...
                '   Pin=', num2str(P.in1(i)), '   Pout=', num2str(P.out(i)),'    Sum=', num2str(P.in1(i)+P.out(i))]);
%                         figure;
%                         plotyy(r_,abs(Ain1),r_,abs(Aout)); title( num2str(round(100*i/NumOfPoints)));
%                         drawnow;
        end
    end
toc    
end

