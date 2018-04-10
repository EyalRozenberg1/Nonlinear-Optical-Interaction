% Sellmeier equation for LBO
% Based on Ghosh, Journal of Applied Physics, 78 (11) 1995

function n = n_lbo_Ghosh(pol, lam, T)

lam_um = lam * 1e6; % switch from meters to micro meters

  if pol=='z'
    A = 1.4489240;
    B = 1.1365228;
    C = 1.1676746e-2;
    D = 1.5830069;
    E = 91;
    G = -446.95031e-6;
%     H = (410.66123 + 1.667e-1*T - 5.1887e-4*T.^2 + 5.5625e-7*T.^3)*1e-6;
    H = 419.33410e-6;
    lamig = 43.5e-9;
  elseif pol == 'x'
    A = 1.4426279;
    B = 1.0109932;
    C = 1.1210197e-2;
    D = 1.2363218;
    E = 91;
    G = -127.70167e-6;
    H = 122.13435e-6;
    lamig = 53e-9;
  elseif pol == 'y'
    A = 1.5014015;
    B = 1.0388217;
    C = 1.2157100e-2;
    D = 1.7567133;
    E = 91;
    G = (372.170 - 2.199e-1*T + 1.1748e-3*T.^2 - 2.05077e-6*T.^3)*1e-6;
%     G = 373.33870e-6;
    H = -415.10435e-6;
    lamig = 32.7e-9;
  else
     error('Bad polarization');
  end


  n     = sqrt(A + B./(1-C./lam_um.^2)+D./(1-E./lam_um.^2));
  R     = lam.^2./(lam.^2-lamig^2);
  dn_dT = (G*R+H*R.^2)./(2*n);
  n     = n + (T-20).*dn_dT;
%   fprintf('n=%f @ lam=%f [um]\n',n(1),lam_um);

end