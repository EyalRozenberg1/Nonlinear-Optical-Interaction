function [A, P] = WavePropagation_SSF2(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples)
% Split step fourier solution for Plane wave propagation

tic 
    switch InteractionType
        case 'Type1 SHG'
            if(PlaneGauss_) % Plane Wave
                P=0;
                % memory allocation
                A.in1 = zeros(NumOfPoints,1);
                A.in2 = 0;
                A.out = zeros(NumOfPoints,1);
                
                % Initialize waves at the start of crystal
                A.in1(1)   = A_from_I(I.in1,n.in1(1));
                A.out(1)   = A_from_I(I.out,n.out(1));
                for i=2:length(CrystalPropAxis)
                    % General Definitions
                    KappaIn  = Kappa(deff,Omega.in1,K.in1,i);
                    Kappaout = Kappa(deff,Omega.out,K.out,i);

                    % coupled equations, boyd pg. 98 equations 2.7.10-12:
                    xi    = i * dx_prop;
                    dAin1 = 2*KappaIn*A.out(i-1).*conj(A.in1(i-1)).*exp(-1i*DeltaK(i)*xi);
                    dAout = Kappaout*A.in1(i-1).^2.*exp(1i*DeltaK(i)*xi);
                    
                    if(Undepleted)
                    A.in1(i) = A.in1(1); % Undepleted Pump
                    else
                    A.in1(i) = A.in1(i-1) + dAin1*dx_prop;
                    end
                    A.out(i) = A.out(i-1) + dAout*dx_prop;
                end
            else % Gaussian Wave
                A         = 0;
                x         = NumOfPoints * dx_prop;
                dx_wide   = dx_prop/2;
                N_wide    = round(x/dx_prop);
                x2        = (-N_wide/2)*dx_wide;
                
                
                % Hankel definitions
                % Note: we can be defined some parameters as const because change in very little
                zr.in1        = pi*n.in1(1)*w0.in1^2/Lambda.in1;     % [m] Reighley range with consider of Lambda in
                w_max.in1     = w0.in1*sqrt(1+((x/2)/zr.in1)^2);
                rmax.max      = 5*w_max.in1;

                mat_H         = hankel_matrix2(0, rmax.max, samples);
                r             = mat_H.r;                  % Radius of beam =	xy coordinate
                fr_sq         = (mat_H.kr/2/pi).^2;       % mat_H.kr	   =	tranform coordinate
                
                % Hankel Transform & Inverse Hankel Transform for out wave
                ht           = @(f) mat_H.JV .* (mat_H.T * (f./mat_H.JR)) /(2*pi);
                iht          = @(F) mat_H.JR .* (mat_H.T * (F./mat_H.JV)) *(2*pi);
                
                % propagation H same as Split Step
                H           = @(PropT, k) exp(PropT)*exp(-1i*k*dx_prop);
                % Prop term
                PropTerm    = @(lambda,n_) 1i*dx_prop*2*pi* realsqrt( (n_/lambda)^2  - fr_sq);
                    
                % according to boyd pg. 118 equations 2.10.5b-c:
                Kin1            = K.in1(1);
                Kout            = K.out(1);
                b_             = @(k,w_0) k*w_0^2;
                b.in1          = b_(Kin1, w0.in1);
                b.out          = b_(Kout, w0.in1);
                
                
                chi_           =  @(B) 2*x2/B;
                chi.in1        = chi_(b.in1);
                chi.out        = chi_(b.out);
                
                % Result for the start of the crystal
                % memory allocation
                P.in1 = zeros(1,NumOfPoints);
                P.out = zeros(1,NumOfPoints);
                
                % according to boyd pg. 117 equation 2.10.5a:
                
                Ain1   = A_from_I(I.in1,n.in1(1))/(1+1i*chi.in1)*exp(-(r.^2)/(w0.in1^2*(1+1i*chi.in1)));
                Aout   = A_from_I(I.out,n.out(1))/(1+1i*chi.out)*exp(-(r.^2)/(w0.in1^2*(1+1i*chi.out)));

                P.in1(1) = P_from_A(Ain1, r, n.in1(1));
                P.out(1) = P_from_A(Aout, r, n.out(1));
                
                for i=2:NumOfPoints
                    
                    % General Definitions
                    KappaIn  = Kappa(deff,Omega.in1,K.in1,i);
                    Kappaout = Kappa(deff,Omega.out,K.out,i);

                    % coupled equations, boyd pg. 98 equations 2.7.10-12:
                    xi    = i * dx_prop;
                    dAin1 = 2*KappaIn*Aout.*conj(Ain1).*exp(-1i*DeltaK(i)*xi);
                    dAout = Kappaout*Ain1.^2.*exp( 1i*DeltaK(i)*xi);
                    
                    Ew_in1 = ht(Ain1);
                    Ew_out = ht(Aout);
                    
                    propterm.in1    = PropTerm(Lambda.in1 ,n.in1(i));
                    propterm.out    = PropTerm(Lambda.out ,n.out(i));
                    
                    Hw.in1 = H(propterm.in1, K.in1(i));
                    Hw.out = H(propterm.out, K.out(i));                    
                
                    Ew_in1 = Hw.in1 .* Ew_in1;
                    Ew_out = Hw.out .* Ew_out;
                    
                    Ain1_temp = iht(Ew_in1);
                    Aout_temp = iht(Ew_out);
                    
                    if(Undepleted)
                        Ain1 = Ain1;
                    else
                        Ain1 = Ain1_temp + dAin1*dx_prop;
                    end
                    Aout = Aout_temp + dAout*dx_prop;
                    
                    if(Undepleted)
                        P.in1(i) = P.in1(1); % Undepleted Pump
                    else
                        P.in1(i) = P_from_A(Ain1, r, n.in1(i));                        
                    end
                    P.out(i) = P_from_A(Aout, r, n.out(i));
                    
                    if(isnan(P.out(i)) || isnan(P.in1(i)))
                        error('ERROR: numerical issue while calculating BW. Exit code');
                    end
                    if(mod(i,round(NumOfPoints/10))==0 || i==NumOfPoints || i==2)
                        disp(['Calc Split step. Done: ', num2str(round(100*i/NumOfPoints)),'%',...
                            '   Pin=', num2str(P.in1(i)), '   Pout=', num2str(P.out(i)),'    Sum=', num2str(P.in1(i)+P.out(i))]);
                        figure;
                        plotyy(r,abs(Ain1),r,abs(Aout)); title( num2str(round(100*i/NumOfPoints)));
                        drawnow;
                    end
                end
            end     
        case 'Type1 THG'
             % compute for this case %
    end
toc    
end

