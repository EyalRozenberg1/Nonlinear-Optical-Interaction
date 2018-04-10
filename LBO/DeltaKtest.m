T_vec   = 140:0.01:160;
Lambda  = 1064.5e-9;

Deltak1 = zeros(size(T_vec));
Deltak2 = zeros(size(T_vec));
Deltak3 = zeros(size(T_vec));

i = 0;
for T = T_vec
    i = i+1;
    Deltak1(i) = DeltaK_CreatorGhosh(T,Lambda);
    Deltak2(i) = 2*2*pi/Lambda*n_lbo_Ghosh('z', Lambda, T) - 2*pi/(Lambda/2)*n_lbo_Ghosh('y', Lambda/2, T);
    Deltak3(i) = 2*2*pi/Lambda*n_lbo_Kato('z', Lambda, T)  - 2*pi/(Lambda/2)*n_lbo_Kato('y', Lambda/2, T);
end

figure; plot(T_vec, Deltak1);
legend('Ghosh1');

figure; plot(T_vec, Deltak1, T_vec, Deltak2, T_vec, Deltak3);
legend('Ghosh1', 'Ghosh2', 'Kato');