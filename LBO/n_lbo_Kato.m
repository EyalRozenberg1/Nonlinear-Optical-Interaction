% Sellmeier equation for LBO
% Based on Kato, IEEE JQE Vol.30 p.2950 (1994).
function n = n_lbo_Kato(pol, lam, T)

 dT     = T-20;
 lam    = lam * 1e6; % switch from meters to micro meters
 
 if strcmp(pol, 'x')
    dn_dT   = (-3.76*lam+2.30)*1e-6;
    dn      = dn_dT.*(dT+29.13e-3*dT.^2);
    n_sq    = 2.4542 + 0.01125./(lam.^2-0.01135) - 0.01388*lam.^2;
 elseif strcmp(pol, 'y')    
    dn_dT   = (6.01*lam-19.40)*1e-6;
    dn      = dn_dT.*(dT-32.89e-4*dT.^2);
    n_sq    = 2.5390 + 0.01277./(lam.^2-0.01189) - 0.01849*lam.^2 + 4.3025e-5*lam.^4 - 2.9131e-5*lam.^6;
 elseif strcmp(pol, 'z')
    dn_dT   = (1.50*lam-9.70)*1e-6;
    dn      = dn_dT.*(dT-74.49*1e-4*dT.^2);
    n_sq    = 2.5865 + 0.01310./(lam.^2-0.01223) - 0.01862*lam.^2 + 4.5778e-5*lam.^4 - 3.2526e-5*lam.^6;
 else
  sprintf('wrong pol (%s)',pol);
 end

    n = sqrt(n_sq) + dn;
end