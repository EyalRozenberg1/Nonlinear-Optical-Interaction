function [A, P] = WavePropagation_SSF(Undepleted, Lambda, w0, NumOfPoints, PlaneGauss_, dx_prop, CrystalPropAxis, DeltaK, K, Omega, n, I, InteractionType, deff, A_from_I, Kappa, P_from_A, samples,NormST)
% Split step fourier solution for Plane wave propagation

tic 
    switch InteractionType
        case 'Type1 SHG'
            if(PlaneGauss_) % Plane Wave
                P=0;
                % memory allocation
                A.in1 = complex(zeros(NumOfPoints,1));
                A.in2 = 0;
                A.out = complex(zeros(NumOfPoints,1));
                
                % Initialize waves at the start of crystal
                A.in1(1)   = A_from_I(I.in1,n.in1(1));
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
                if NormST==0
                    A=0;
                end
                x_max       = NumOfPoints * dx_prop;
                % Hankel definitions
                % Note: we can be defined some parameters as const because change in very little
                Kin1        = K.in1(1);
                
                zr.in1      = pi*n.in1(1)*w0.in1^2/Lambda.in1;     % [m] Reighley range with consider of Lambda in
                w_max.in1   = w0.in1*sqrt(1+(x_max/zr.in1)^2);
                
                rmax.in1      = 2*w_max.in1;
                rmax.max      = rmax.in1;

                mat_H2        = hankel_matrix2(0, rmax.max, samples);
                r_            = mat_H2.r;                  % Radial co-ordinate vector
                fr_sq_        = (mat_H2.v).^2;             % frequency co-ordinate vector
                
                % Hankel Transform & Inverse Hankel Transform for out wave
                ht2         = @(f) mat_H2.JV .* (mat_H2.T * (f./mat_H2.JR)) /(2*pi);
                iht2        = @(F) mat_H2.JR .* (mat_H2.T * (F./mat_H2.JV)) *(2*pi);
                
                % propagation H same as Split Step
                H              = @(PropT, k)  exp(PropT).*exp(-1i*k*dx_prop);
                PropTerm       = @(lambda,n_) 1i*dx_prop*2*pi*sqrt((n_/lambda)^2 - fr_sq_);
                
                % according to boyd pg. 118 equations 2.10.5b-c:
                b_             = @(k,w_0) k*w_0^2;
                b.in1          = b_(Kin1, w0.in1);
                
                x_waist        = x_max/2; % [m] waists location in the crystal
                chi_           = @(B) 2*(0-x_waist)/B;
                chi.in1        = chi_(b.in1);
                
                % Result for the start of the crystal
                % memory allocation
                P.in1 = complex(zeros(1,NumOfPoints));
                P.out = complex(zeros(1,NumOfPoints));
                
                % according to boyd pg. 117 equation 2.10.5a:
                Ain1        = A_from_I(I.in1,n.in1(1))/(1+1i*chi.in1)*exp(-(r_.^2)/(w0.in1^2*(1+1i*chi.in1)));
                Aout        = complex(zeros(samples,1));
                
%                 dAin1Phase   = zeros(NumOfPoints,1);
%                 dAoutPhase   = zeros(NumOfPoints,1);
%                 Ain1PhaseTemp   = zeros(NumOfPoints,1);
%                 AoutPhaseTemp   = zeros(NumOfPoints,1);

%                 Ain1Phase   = zeros(NumOfPoints,1);
%                 AoutPhase   = zeros(NumOfPoints,1);

                P.in1(1) = P_from_A(Ain1, r_, n.in1(1));
                P.out(1) = 0;
                
                for i=1:NumOfPoints
                    
                    % Prop term
                    propterm.in1    = PropTerm(Lambda.in1 ,n.in1(i));
                    propterm.out    = PropTerm(Lambda.out ,n.out(i));
                    
                    propterm.in1( real(propterm.in1)>0 ) = -propterm.in1( real(propterm.in1)>0 );
                    propterm.out( real(propterm.out)>0 ) = -propterm.out( real(propterm.out)>0 );

                    Ew_in1 = ht2(Ain1);
                    Ew_out = ht2(Aout);
                    
                    Hw.in1 = H(propterm.in1, K.in1(i));
                    Hw.out = H(propterm.out, K.out(i));                    
                
                    Ew_in1 = Ew_in1 .* Hw.in1;
                    Ew_out = Ew_out .* Hw.out;
                    
                    Ain1_temp = iht2(Ew_in1);
                    Aout_temp = iht2(Ew_out);
                    
                    % coupled equations, boyd pg. 98 equations 2.7.10-12:
                    % General Definitions
                    KappaIn  = Kappa(deff,Omega.in1,K.in1,i);
                    Kappaout = Kappa(deff,Omega.out,K.out,i);
                    
                    xi    = i * dx_prop;% changed
                    dAin1 = 2*KappaIn*Aout.*conj(Ain1)*exp( -1i*DeltaK(i)*xi );
                    dAout = Kappaout*Ain1.^2*exp( 1i*DeltaK(i)*xi );
                    
                    if(Undepleted)
                        Ain1 = Ain1;
                    else
                        Ain1 = Ain1_temp + dAin1*dx_prop;
                    end
                    Aout = Aout_temp + dAout*dx_prop;

%                     Ain1PhaseTemp(i) = unwrap(angle(Ain1_temp(1)));
%                     AoutPhaseTemp(i) = unwrap(angle(Aout_temp(1)));
%                     dAin1Phase(i) = unwrap(angle(dAin1(1)));
%                     dAoutPhase(i) = unwrap(angle(dAout(1)));
                    
%                     Ain1Phase(i) = unwrap(angle(Ain1(1)));
%                     AoutPhase(i) = unwrap(angle(Aout(1)));
                    
                    if(Undepleted)
                        P.in1(i) = P.in1(1); % Undepleted Pump
                    else
                        P.in1(i) = P_from_A(Ain1, r_, n.in1(i));                        
                    end
                    P.out(i) = P_from_A(Aout, r_, n.out(i));
                    
                    if(isnan(P.out(i)) || isnan(P.in1(i)))
                        error('ERROR: numerical issue while calculating BW. Exit code');
                    end
                    if(mod(i,round(NumOfPoints/10))==0 || i==NumOfPoints || i==2)
                        disp(['Calc Split step. Done: ', num2str(round(100*i/NumOfPoints)),'%',...
                            '   Pin=', num2str(P.in1(i)), '   Pout=', num2str(P.out(i)),'    Sum=', num2str(P.in1(i)+P.out(i))]);
%                         figure;
% %                         plotyy(r_,abs(Ain1),r_,abs(Aout)); title( num2str(round(100*i/NumOfPoints)));
%                         title( ['\zeta=',num2str(round(100*i/NumOfPoints)),'%L']);
%                         yyaxis left
% %                         plot(r_,abs(Aout)./max(abs(Aout)),'b','linewidth',2);
%                         plot([-r_ r_],[abs(Aout)./max(abs(Aout))  abs(Aout)./max(abs(Aout))],'b-','linewidth',2);
%                         xlabel('transverse axis [m]')
%                         ylabel('normalized |A_{2\omega}|')
%                         yyaxis right
% %                         plot(r_,abs(Ain1)./max(abs(Ain1)),'r','linewidth',2);
%                         plot([-r_ r_],[abs(Ain1)./max(abs(Ain1))  abs(Ain1)./max(abs(Ain1))],'r-','linewidth',2);
%                         ylabel('normalized |A_{\omega}|')
%                         drawnow;
                    end
                end
                if NormST==1
                    A.in1 = A_from_I(P.in1.'/(pi*w0.I^2),n.in1);
                    A.in2 = 0;
                    A.out = A_from_I(P.out.'/(pi*w0.I^2),n.out);
                end
            end     
        case 'Type1 THG'
             % compute for this case %
    end
toc
% figure; 
% % plot((Ain1Phase));
% % hold on;
% % plot((AoutPhase));
% plot((2*AoutPhase-Ain1Phase)-mean((2*AoutPhase-Ain1Phase)));
% hold on;
% plot(-atan((CrystalPropAxis'-x_max/2)/zr.in1))
% % legend('Guoy(A\omega)','Guoy(A2\omega)','Guoy(A2\omega)-Guoy(A\omega)')
% legend('2*Guoy(A2\omega)-Guoy(A\omega)','-atan(z/z_r)')
% keyboard;
end

