function [A, P] = WavePropagation_SSF_first_version(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples)
% Split step fourier solution for Plane wave propagation
tic
    A=0;
%     P=0;
    switch InteractionType
        case 'Type1 SHG'
            if(PlaneGauss_) % Plane Wave
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
                    dAout = Kappaout*A.in1(i-1).^2.*exp( 1i*DeltaK(i)*xi);
                    
                    if(Undepleted)
                    A.in1(i) = A.in1(1); % Undepleted Pump
                    else
                    A.in1(i) = A.in1(i-1) + dAin1*dx_prop;
                    end
                    A.out(i) = A.out(i-1) + dAout*dx_prop;
                end
            else % Gaussian Wave
                % Hankel definitions
                % Note: we can be defined some parameters as const because change in very little
                nConst      = n.in1(1);
                Kin1        = K.in1(1);
                Kout        = K.out(1);
                
                zr.in1      = pi*n.in1(1)*w0.in1^2/Lambda.in1;     % Reighley range with consider of Lambda in
                %zr.out      = zr.in1;
                zr.out      = pi*n.out(1)*w0.out^2/Lambda.out;     % Reighley range with consider of Lambda out
                
                w_max.in1   = w0.in1*sqrt(1+(NumOfPoints*dx_prop/zr.in1)^2);
                %w_max.out   = w_max.in1;
                w_max.out   = w0.out*sqrt(1+(NumOfPoints*dx_prop/zr.out)^2);
                
                rmax.in1        = 2*w_max.in1;
                rmax.out        = 2*w_max.out;
                rmax.max        = 1.2 * max(rmax.in1, rmax.out);

                mat_H2_in1      = hankel_matrix2(0, rmax.max, samples);
                r_in1           = mat_H2_in1.r;                  % Radius of beam =	xy coordinate
                fr_sq_in1       = (mat_H2_in1.v).^2;             %   =	tranform coordinate
                
                % Hankel Transform & Inverse Hankel Transform for in wave
                ht2_in1         = @(f) mat_H2_in1.JV .* (mat_H2_in1.T * (f./mat_H2_in1.JR)) /(2*pi);
                iht2_in1        = @(F) mat_H2_in1.JR .* (mat_H2_in1.T * (F./mat_H2_in1.JV)) *(2*pi);
                
                %mat_H2_out      = mat_H2_in1;
                mat_H2_out      = hankel_matrix2(0, rmax.max, samples);
                %r_out           = r_in1;
                r_out           = mat_H2_out.r;                  % Radius of beam =	xy coordinate
                %fr_sq_out       = fr_sq_in1;
                fr_sq_out       = (mat_H2_out.v).^2;             %   =	tranform coordinate
                
                % Hankel Transform & Inverse Hankel Transform for out wave
                %ht2_out         = ht2_in1;
                ht2_out         = @(f) mat_H2_out.JV .* (mat_H2_out.T * (f./mat_H2_out.JR)) /(2*pi);
                %iht2_out        = iht2_in1;
                iht2_out        = @(F) mat_H2_out.JR .* (mat_H2_out.T * (F./mat_H2_out.JV)) *(2*pi);
                
                % propagation H same as Split Step
                % Transmition Function
                H              = @(PropT, k) exp(PropT).*exp(-1i*k*dx_prop);
                
                b_             = @(k,w_0) k*w_0^2;
                b.in1          = b_(Kin1, w0.in1);
                %b.out          = b.in1;
                b.out          = b_(Kout, w0.out);
                
                chi_           = @(B) 2*((0-NumOfPoints/2)*dx_prop)./B;
                chi.in1        = chi_(b.in1);
                chi.out        = chi_(b.out);
                
                % Result for the start of the crystal
                % memory allocation
                P.in1 = zeros(1,NumOfPoints);
                P.out = zeros(1,NumOfPoints);
                
                Ain1   = A_from_I(I.in1,n.in1(1))./(1+1i*chi.in1).*exp(-r_in1.^2/w0.in1^2./(1+1i*chi.in1));
                Aout   = A_from_I(I.out,n.out(1))./(1+1i*chi.out).*exp(-r_out.^2/w0.out^2./(1+1i*chi.out));

                P.in1(1) = P_from_A(Ain1, r_in1, n.in1(1));
                P.out(1) = P_from_A(Aout, r_out, n.out(1));
                %f=figure;
                
                for i=2:NumOfPoints
                    % Prop term
                    PropTerm.in1    = @(lambda,fr_sq) 1i*dx_prop*2*pi* sqrt( ((n.in1(i)/lambda).^2 - fr_sq));
                    propterm.in1    = PropTerm.in1(Lambda.in1, fr_sq_in1);
                    propterm.in1( real(propterm.in1)>0 ) = -propterm.in1( real(propterm.in1)>0 ) ;
                    
                    PropTerm.out    = @(lambda,fr_sq) 1i*dx_prop*2*pi* sqrt( ((n.out(i)/lambda).^2 - fr_sq));
                    propterm.out    = PropTerm.out(Lambda.out, fr_sq_out);
                    propterm.out( real(propterm.out)>0 ) = -propterm.out( real(propterm.out)>0 ) ;
                    % General Definitions
                    
                    KappaIn  = Kappa(deff,Omega.in1,K.in1,i);
                    Kappaout = Kappa(deff,Omega.out,K.out,i);

                    % coupled equations, boyd pg. 98 equations 2.7.10-12:
                    xi    = i * dx_prop;
                    dAin1 = 2*KappaIn*Aout.*conj(Ain1).*exp(-1i*DeltaK(i)*xi);
                    dAout = Kappaout*Ain1.^2.*exp( 1i*DeltaK(i)*xi);
                    
                    Ew_in1 = ht2_in1(Ain1);
                    Ew_out = ht2_out(Aout);
                    
                    %Hw.in1 = H(propterm.in1, K.in1(1));
                    Hw.in1 = H(propterm.in1, K.in1(i));
                    %Hw.out = H(propterm.out, K.out(1));
                    Hw.out = H(propterm.out, K.out(i));                    
                
                    Ew_in1 = Ew_in1 .* Hw.in1;
                    Ew_out = Ew_out .* Hw.out;
                    
                    Ain1_temp = iht2_in1(Ew_in1);
                    Aout_temp = iht2_out(Ew_out);
                    
                    if(Undepleted)
                        Ain1 = Ain1;
                    else
                        Ain1 = Ain1_temp + dAin1*dx_prop;
                    end
                    Aout = Aout_temp + dAout*dx_prop;
                    
                    %P.in1(i) = P_from_A(Ain1, r_in1, n.in1(1));
                    P.in1(i) = P_from_A(Ain1, r_in1, n.in1(i));
                    
                    %P.out(i) = P_from_A(Aout, r_out, n.out(1));
                    P.out(i) = P_from_A(Aout, r_out, n.out(i));
                    
                    if(isnan(P.out(i)) || isnan(P.in1(i)))
                        error('ERROR: numerical issue while calculating BW. Exit code');
                    end
                    if(mod(i,round(NumOfPoints/10))==0 || i==NumOfPoints || i==2)
                        disp(['Calc Split step. Done: ', num2str(round(100*i/NumOfPoints)),'%',...
                            '   Pin=', num2str(P.in1(i)), '   Pout=', num2str(P.out(i)),'    Sum=', num2str(P.in1(i)+P.out(i))]);
%                         plotyy(r_in1,abs(Ain1),r_in1,abs(Aout)); title( num2str(round(100*i/NumOfPoints)));
%                         drawnow;
                    end
                end
            end     
        case 'Type1 THG'
             % compute for this case %
    end
toc    
end

