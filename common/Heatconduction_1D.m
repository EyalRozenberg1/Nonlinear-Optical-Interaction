%%% 1D Heat conduction
%%% Martin Albers
%%%[1] Strang G. and Fix G. (2008): An analysis of the Finite Element Method,Second Edition, Wellesley-Cambridge Press, Wellesley USA
%%%[2] R. W. Lewis et al. (1996):  The Finite Element Method in Heat Transfer Analysis, John Wiley and Sons, West Sussex England

%%%--- Problem ---%%%
% Solving 1-D Heat diffusion in a unit rod
% 
clear all
clc

%%%--- Mesh generation ---%%%
h       = 1e-4;                   % Distance between nodes [m]
units   = 5e-2;                   % Length of rod
p       = [0:h:units]';           % Nr of nodes
x       = [0:h:units];            % Spatial discretisation
nr_nodes = length(p);             % Nr of nodes

  for i = 2:nr_nodes        % Creates the element matrix
      e(i-1,:) = [i-1 i];
  end
  
%%%--- Creating space ---%%%
nr_elements =length(e);     % Nr of elements

K=zeros(nr_nodes,nr_nodes); % Empty stiffness matrix
M=zeros(nr_nodes,nr_nodes); % Empty mass matrix
F=zeros(nr_nodes,1);        % Empty forcing vector
Q=zeros(nr_nodes,1);        % Initial forcing vector, 

%%%--- Parameters ---%%%
P=1;                       % Heat conductivity
pc=1;                      % Heat capacity
A=1;                       % Element cross section

%%%--- Assembling matrices ---%%%
for k = 1:nr_elements                        % Loop over all elements
    
        nodes=e(k,:);                        % Identify nodes
        h=p(k+1)-p(k);                       % Calculate distance between nodes
        
        Ke = ((A*P)/h)*[1 -1;-1 1];          % Stiffness entry
        Me = ((pc*A*h)/6)*[2 1;1 2];         % Mass entry       
        Fe = ((A*h)/6)*[2 1;1 2]*Q(nodes,1); % Forcing entry
                            
        K(nodes,nodes)=K(nodes,nodes)+Ke;    % Position stiffness entry at identified nodes
        M(nodes,nodes)=M(nodes,nodes)+Me;    % Position mass entry at identified nodes      
        F(nodes,1)=F(nodes,1)+Fe;            % Position forcing entry at identified nodes          
end


%%%--- Initial conditions ---%%%
Told=zeros(nr_nodes,1);     % Initial temperature, 0 everywhere
F(1)=1;                     % With unit heat production at node one

%%%--- Time parameters ---%%%
dt = 0.01;                  % Size of time step
t=1;                        % Time

%%%--- Preparing matrices ---%%%
% See [2] page 108

A=((1/2)*K+(1/(dt))*M);     
B=(-(1/2)*K+M/(dt));        

C=(M+0.5*dt*K);             
D=(M-0.5*dt*K);             

R=((1/3)*K+(1/(2*dt))*M);   
S=((-1/6)*K+(M/(2*dt)));    
 
R(nr_nodes,:)=0;
R(:,nr_nodes)=0;
R(nr_nodes,nr_nodes)=1;

S(nr_nodes,:)=0;
S(:,nr_nodes)=0;
S(nr_nodes,nr_nodes)=1;

Tnew = Told;
Fold = F;
Fnew = Fold;

%%%--- Time stepping ---%%%
for i = 1:t/dt
    
    thermometer(i)=Tnew(1);                    % Vector to store temperature data at node 1

  %  Tnew = C\(D*Told+dt*F) ;                  % Crank Nicolson
   Tnew=R\((S*Told+((1/6)*(Fold+2*Fnew))));    % Galerking
 
   Told=Tnew;

end


%%%--- Analytical solution ---%%%
tt=[0:dt:t-dt];         % Time discretisation
        
for nt = 1:length(tt)   % Loop over time
 
% Analytical solution
Tana = 2*(tt(nt)/pi)^(1/2)*((exp((-x.^2./(4*tt(nt))))-0.5.*x.*(sqrt(pi/tt(nt)).*erfc((x./(2*sqrt(tt(nt))))))));

thermometerTana(nt)=Tana(1);    % Store temperature at node 1
end

  figure(1) % [2] page 109
  plot(x,Tnew,'ro',x,Tana,'b')
  title('Temperature distribution in rod')
  xlabel('Length')
  ylabel('Temperature [C]')
  legend(['Galerkin dt=',num2str(dt)],'Exact')
  axis([0 2 0 1.4])
%     
  figure(2) % [2] page 110
  plot(tt,thermometerTana,tt,thermometer,'r')
  title('Temperature at point x=0 (Heating Curve)')
  xlabel('Time [s]')
  ylabel('Temperature [C]')
  legend('Exact',['Galerkin dt=', num2str(dt)])
%     
    
%  Tnew1=Tnew;
%  Tana1=Tana;

%   Tnew05=Tnew;
%   Tana05=Tana;
% %     
%   Tnew01=Tnew;
%   Tana01=Tana;
%   
%   
%   figure(1)
%   plot(x,Tnew,'ro',x,Tana,'b',x,Tnew1,'ro',x,Tana1,'b',x,Tnew05,'ro',x,Tana05,'b')
%   title('Temperature distribution in rod')
%   xlabel('Length')
%   ylabel('Temperature [C]')
%   legend('Galerkin dt=0.01','Exact')
%   axis([0 2 0 1.4])