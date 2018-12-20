function [figNum] = sim_and_stationary_states_norm_units_for_JOSAB2(CrystalPropAxis, A, Omega, K, DeltaK, deff, c, figNum)

font_name = 'times';
line_width = 2;
q1_color = [ 1 0 0 ];
q3_color = [ 0 0 1 ];


gamma   = @(omg,k) 2*deff*omg^2./(k*c^2); % gamma = @(omg,k) deff*omg^2./(k*c^2);
gamma1  = gamma(Omega.in1,K.in1);
gamma2  = gamma(Omega.in2,K.in2);
gamma3  = gamma(Omega.out,K.out);

A1 = A.in1(:,1);
A2 = A.in1(:,1);
A3 = A.out(:,1);

q1_sq = abs(A1./sqrt(gamma1)).^2;
q2_sq = abs(A2./sqrt(gamma2)).^2;
q3_sq = abs(A3./sqrt(gamma3)).^2;

DeltaGamma = DeltaK./sqrt(gamma1.*gamma2.*gamma3);
xi = CrystalPropAxis'.*sqrt(gamma1.*gamma2.*gamma3);

K1 = q1_sq + q3_sq;
% K2 = q1_sq - q2_sq;
K3 = q2_sq + q3_sq;

P1 = q1_sq + q2_sq - 2*q3_sq;
P3 = K1+K3;
theta_minus = (5*DeltaGamma-sqrt(DeltaGamma.^2+6*P3))/6;

q1_minus = sqrt((DeltaGamma-theta_minus).*(DeltaGamma-2*theta_minus)).*exp(1i*theta_minus.*xi);
q1_minus(DeltaGamma>sqrt(2*P3))=0;

% q2_minus = q1_minus;

q3_minus = (DeltaGamma-theta_minus).*exp(2*1i*theta_minus.*xi);
q3_minus(DeltaGamma>sqrt(2*P3)) = sqrt(P3(DeltaGamma>sqrt(2*P3))/2).*exp(1i*DeltaGamma(DeltaGamma>sqrt(2*P3)).*xi(DeltaGamma>sqrt(2*P3)));

% q stats
% --------------------------------------------------------- %
figure;
hold on;
set(gcf,'color','white');
plot(DeltaGamma./sqrt(P3), q1_sq./P3, '.', 'color', q1_color );
plot(DeltaGamma./sqrt(P3), q3_sq./P3, '.', 'color', q3_color );
plot(DeltaGamma./sqrt(P3), abs(q1_minus).^2./P3, 'k--', 'color', q1_color );
plot(DeltaGamma./sqrt(P3), abs(q3_minus).^2./P3, 'k--', 'color', q3_color );
hold off;
set(gca, 'FontName', font_name);
xlabel('\Delta\Gamma/P_3^{1/2}');
ylabel('|q|^2/P_3');
legend( '|q_1|^2', '|q_3|^-', 'q_1^-', 'q_3^-');

% --------------------------------------------------------- %


P1_minus = 2*abs((DeltaGamma-theta_minus).*(DeltaGamma-2*theta_minus))-2*(DeltaGamma-theta_minus).^2;
% calculate and plot P1
% --------------------------------------------------------- %
figure;
hold on;
set(gcf,'color','white');
plot(DeltaGamma./sqrt(P3), P1./P3, '.', 'color', q3_color);
plot(DeltaGamma./sqrt(P3), P1_minus./P3, 'k--', 'LineWidth', line_width);
hold off;
set(gca, 'FontName', font_name);
ylabel('P_1/P_3');
xlabel('\Delta\Gamma/P_3^{1/2}');
legend('P_1','P_-');


% calculate rnl 1
% --------------------------------------------------------- %
d_P1 = P1-P1_minus;
figure;
hold on;
set(gcf,'color','white');
plot(xi.*sqrt(P3), abs(d_P1./P3), '.', 'color', q3_color );
% plot([mean(xi.*sqrt(P3)) mean(xi.*sqrt(P3))], [0 1], 'k--' );
hold off;
set(gca, 'FontName', font_name);
xlabel('\xiP_3^{1/2}');
ylabel('r_{nl}1');


% calculate rnl 2
% --------------------------------------------------------- %
figure;
hold on;
set(gcf,'color','white');
rnl_2 = 2/sqrt(27)*1./P3(2:end).*abs((DeltaGamma(2:end)-DeltaGamma(1:end-1))./(xi(2:end)-xi(1:end-1)));
plot(xi(2:end).*sqrt(P3(2:end)), rnl_2, '.', 'color', q3_color );
ylim([0 max(rnl_2)])
% plot([mean(xi.*sqrt(P3)) mean(xi.*sqrt(P3))], [0 max(rnl_2)], 'k--' );
hold off;
set(gca, 'FontName', font_name);
xlabel('\xiP_3^{1/2}');
ylabel('r_{nl}')

% calculate rnl 3 - not normalized
% --------------------------------------------------------- %
% figure;
% hold on;
% set(gcf,'color','white');
% plot(CrystalPropAxis(2:end), rnl_2, '.', 'color', q3_color );
% % plot([mean(CrystalPropAxis) mean(CrystalPropAxis)], [0 max(rnl_2)], 'k--' );
% hold off;
% set(gca, 'FontName', font_name);
% xlabel('\zeta[m]');
% ylabel('r_{nl}2')
