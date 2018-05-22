
N   = 50; % Number of grid in x,y-direction
L1  = 5e-2; % Domain size
L2  = 5e-3; % Domain size
% Grid point
x = linspace(0,L1,N);
y = linspace(0,L2,N);
% Make it staggered.
x = (x(1:end-1)+x(2:end))/2;
y = (y(1:end-1)+y(2:end))/2;
[X,Y] = meshgrid(x,y);

% Let's use MATLAB logo.
% A variable u0 is defined at the center of each grid cell
% thus the number of grid point is N-1.
u0(:,:) = 273.15 + 150*ones(49);%peaks(N-1);
% Plot it
handle_surf = surf(X,Y,u0);
handle_axes = gca;
handle_axes.ZLim = 273.15 + [120,180];
handle_axes.CLim = 273.15 + [120,180];
title('Evolution of MATLAB Logo by Heat equation');

dx = x(2)-x(1); % spatial grid size
alpha = 10.8e-5; % [1/K] coefficient
tspan = linspace(130,170,40);
[t,u] = ode15s(@(t,x)getRHS(x,alpha,dx,N),tspan,u0(:));

Tn = length(t);
u = reshape(u,Tn,N-1,N-1);
filename = 'heat.gif';
for ii=1:Tn
    Z = u(ii,:,:);
    Z = squeeze(Z);
    handle_surf.ZData = Z;
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if ii==1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end
