
% article Gil Porat Fully Nonlinear Adiabatic JOSAB (2013)- Gil
% "Efficient, broadband, and robust frequency conversion by fully nonlinear adiabatic three-wave mixing"

function [figNum] = sim_and_stationary_states_norm_units_for_JOSAB(CrystalPropAxis, A, Omega, K, DeltaK, deff, c, figNum)

font_name = 'times';
font_size = 18;
line_width = 2;
q1_color = [ 1 0 0 ];
q3_color = [ 0 0 1 ];

gamma   = @(omg,k) 2*deff*omg^2./(k*c^2); % gamma = @(omg,k) deff*omg^2./(k*c^2);
gamma1  = gamma(Omega.in1,K.in1);
gamma2  = gamma(Omega.in2,K.in2);
gamma3  = gamma(Omega.out,K.out);

q           = @(a,gmma) a./sqrt(gmma); % generalized coordinates
q1          = q(A.in1(:,1).',gamma1);
q2          = q1;
q3          = q(A.out(:,1).',gamma3);

q1_sq       = abs(q1).^2;
q2_sq       = abs(q2).^2;
q3_sq       = abs(q3).^2;

DeltaGamma  = DeltaK./sqrt(gamma1.*gamma2.*gamma3); % the relative strength of the phase mismatch compared to the nonlinearity
xi          = CrystalPropAxis.*sqrt( gamma1.*gamma2.*gamma3); % the scaled propagation length

K1          = q1_sq + q3_sq;
K2          = q1_sq - q2_sq;
K3          = q2_sq + q3_sq;

P1          = q1_sq + q2_sq - 2*q3_sq;
P2          = K2;% P2 equal to zero while the two lower frequency has the same photon flux
P3          = K1 + K3;

% plot the minus solution
theta_sign  = -1;
theta       = 1/6 * ( 5*DeltaGamma + theta_sign * sqrt( DeltaGamma.^2 + 6*P3 ) );

inds_below = find( DeltaGamma < sqrt(2*P3) );

q1_minus             = zeros( 1, length( DeltaGamma ) );
% q1_minus(inds_below) = sqrt( (DeltaGamma(inds_below) - theta(inds_below)).*(DeltaGamma(inds_below) - 2*theta(inds_below)) );  Changed at 2017_06_20
q1_minus(inds_below) = sqrt( (DeltaGamma(inds_below) - theta(inds_below)).*(DeltaGamma(inds_below) - 2*theta(inds_below)) ).* exp( 1i*theta(inds_below).*xi(inds_below) );

q2_minus = q1_minus;

% q3_minus             = sqrt(P3/2) .* ones( 1, length( DeltaGamma )); Changed at 2017_06_20
% q3_minus(inds_below) = DeltaGamma(inds_below) - theta(inds_below);
q3_minus             = sqrt(P3/2) .* exp(1i*DeltaGamma.*xi);
q3_minus(inds_below) = ( DeltaGamma(inds_below) - theta(inds_below) ).* exp( 2i*theta(inds_below).*xi(inds_below) );

figNum = figNum + 1;
figure(figNum);set(gcf,'color','white');
hold on
plot( DeltaGamma./sqrt(P3), q1_sq./P3, '.', 'color', q1_color );
plot( DeltaGamma./sqrt(P3), q3_sq./P3, '.', 'color', q3_color );
set(gca, 'FontName', font_name);
axis( [ min( DeltaGamma./sqrt(P3) ) max( DeltaGamma./sqrt(P3) ) 0 max(q1_sq./P3)*1.25 ] );

% expected results 
plot( DeltaGamma./sqrt(P3), (abs(q1_minus).^2)./P3, 'k--', 'color', q1_color );
plot( DeltaGamma./sqrt(P3), (abs(q3_minus).^2)./P3, 'k--', 'color', q3_color );
set(gca, 'FontName', font_name);
xlabel('\Delta\Gamma/P_3^{1/2}');
ylabel('|q|^2/P_3');
legend( '|q_1|^2', '|q_3|^-', 'q_1^-', 'q_3^-', 'Location', 'NorthWest' );


% calculate r, the adiabatic criterion parameter
% --------------------------------------------------------- % got to here with readding at 20/06/2017
% DG_min  = -10*P3;
% DG_max  =  10*P3;
% L_xi    = 1./sqrt(P3);
% 
% d_theta_d_Gamma = 1/6 * ( 5 + theta_sign * DeltaGamma ./ sqrt( DeltaGamma.^2 + 6*P3 ) );
% 
% sign1           = sign( ( DeltaGamma - 2.*theta ).* ( DeltaGamma - theta ) );
% dP1_0_d_Gamma   = sign1 .* 2 .* ( ( 1 - d_theta_d_Gamma ) .* ( DeltaGamma - 2*theta ) + ( 1 - 2*d_theta_d_Gamma ) .* ( DeltaGamma - theta ) ) - 2 .* ( 2*( DeltaGamma - theta ) .* ( 1 - d_theta_d_Gamma ) );
% 
% d2H_dQ1         = 64 * abs( ( DeltaGamma - theta ).^2 .* ( DeltaGamma - 2*theta ) );
% 
% d2H_dP1         =   -1 ./ ( 8 * abs( DeltaGamma - theta ) ) + ...
%                     1 ./ ( 16 * abs( DeltaGamma - 2*theta ) ) - ...
%                     DeltaGamma ./ ( 2 * abs( (DeltaGamma - theta ) .* (DeltaGamma - 2*theta ) ) ) + ...
%                     DeltaGamma ./ ( 4 * abs( DeltaGamma - theta ).^2 );
% 
% nu              = 1i * sqrt( d2H_dQ1 .* d2H_dP1 );
% d_DG_dxi        = ( DG_max - DG_min ) ./ L_xi; % TODO: change to nonlinear as well
% r               = abs(1./nu .* dP1_0_d_Gamma .* d_DG_dxi./ P3); %1./abs(nu) .* abs( dP1_0_d_Gamma ) .* d_DG_dxi ./ P3

r = abs(...
        ((0.5*(abs(q1_sq)+abs(q2_sq))-abs(q3_sq)) - ...
         (0.5*(abs(q1_minus).^2+abs(q2_minus).^2) - abs(q3_minus).^2))./...
        (0.5*(abs(q1_sq)+abs(q2_sq))+abs(q3_sq))...
        );
figNum = figNum + 1;
figure(figNum);set(gcf,'color','white');
plot( xi.*sqrt(P3), r);
set(gca, 'FontName', font_name);
xlabel('\xi');
ylabel('r_{nl}');
FormatPlotFontSizeNameLine( gca, gcf, font_size, 'times', line_width );
axis( [ min( xi.*sqrt(P3) ) max( xi.*sqrt(P3) ) 0 max(r)*1.05 ] );
[M,I] = max(r);
xi_mul_P3 = xi.*sqrt(P3);
text(xi_mul_P3(I), r(I), ['\leftarrow ', num2str(M)]);

% calculate and plot P1
P1_minus = 2*abs( (DeltaGamma - theta).*(DeltaGamma - 2*theta) ) - 2*(DeltaGamma - theta).^2;
figNum = figNum + 1;
figure(figNum)
hold on;
% plot( xi.*sqrt(P3), P1./P3); Changed at 2017_06_20
plot( DeltaGamma.*sqrt(P3), P1./P3);
% plot( xi.*sqrt(P3), P1_minus./P3, 'k--', 'LineWidth', line_width ); Changed at 2017_06_20
plot( DeltaGamma.*sqrt(P3), P1_minus./P3, 'k--', 'LineWidth', line_width );
set(gca, 'FontName', font_name);
FormatPlotFontSizeNameLine( gca, gcf, font_size, 'times', line_width );
ylabel('P_1/P_3');
% xlabel('\xi'); Changed at 2017_06_20
xlabel('\Delta\Gamma/P_3^{1/2}');
legend('P_1','P_-');
hold off;
