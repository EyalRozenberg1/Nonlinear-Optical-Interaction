% Load saved figures
a=hgload('eta_vs_Pw@137.fig');
b=hgload('eta_vs_Pw@141.fig');
c=hgload('eta_vs_Pw@145.fig');
d=hgload('eta_vs_Pw@149(vs 147 simulation).fig');
% Prepare subplots
figure
h(1)=subplot(2,2,1);
h(2)=subplot(2,2,2);
h(3)=subplot(2,2,3);
h(4)=subplot(2,2,4);
% Paste figures on the subplots
copyobj(allchild(get(a,'CurrentAxes')),h(1));
copyobj(allchild(get(b,'CurrentAxes')),h(2));
copyobj(allchild(get(c,'CurrentAxes')),h(3));
copyobj(allchild(get(d,'CurrentAxes')),h(4));
% Add legends
l(1)=title(h(1),'137[^\circC] ');
l(2)=title(h(2),'141[^\circC]');
l(3)=title(h(3),'145[^\circC]');
l(4)=title(h(4),'149[^\circC]');

ylabel(h(1),'\eta');
ylabel(h(2),'\eta');
ylabel(h(3),'\eta');
ylabel(h(4),'\eta');

xlabel(h(1),'P\omega[W]');
xlabel(h(2),'P\omega[W]');
xlabel(h(3),'P\omega[W]');
xlabel(h(4),'P\omega[W]');