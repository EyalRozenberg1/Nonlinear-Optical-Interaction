xi              = 5.1319;
L               = 5e-2;
zr              = L/(2*xi);
dx_prop         = 5e-5;

MatDim      = 200;

Lambda.in1  = 1064.2e-9;   % Wavelength [m]
Lambda.in2  = Lambda.in1;
Lambda.out  = Lambda.in1/2;

Pol.in1     = 'z';
Pol.in2     = 'z';
Pol.out     = 'y';

c           = 299792458;               % m/s Speed of light
eps0        = 8.854187817e-12;         % F/m
k           = @(n, lambda)              n.*2*pi/lambda;
w           = @(lambda)                 2*pi*c/lambda;
refIdx      = @(pol, lambda, temp)      n_lbo_Kato(pol, lambda, temp);

NumOfPoints     = round(L/dx_prop);
% Z = zeros(MatDim,MatDim,NumOfPoints);
Z = linspace(0,L,NumOfPoints);
Z = Z + Z(2);

LL = zeros(MatDim,MatDim,NumOfPoints);
ll = ( (atan((L/2)/zr) + atan((Z-L/2)/zr))./Z - mean((atan((L/2)/zr) + atan((Z-L/2)/zr))./Z) ) ...
/ max(atan((L/2)/zr) + atan((Z-L/2)/zr)./Z);

for i=1:MatDim
    for j=1:MatDim
        LL(i,j,:) = ll;
    end
end


Tpm = repmat(linspace(148,150,MatDim)',1,MatDim,NumOfPoints);
Dt  = repmat(linspace(-1,1,MatDim),MatDim,1,NumOfPoints);

TempGrad = Tpm + 0.5*(Dt) .* LL;

GouyPhase = -atan((Z-L/2)/zr).';
DK = zeros(MatDim,MatDim,NumOfPoints);
obj = zeros(MatDim,MatDim,NumOfPoints);
std_ = zeros(MatDim,MatDim);
mean_ = zeros(MatDim,MatDim);

for i=1:MatDim
    for j=1:MatDim
        DK(i,j,:)  = DeltaK_Creator(TempGrad(i,j,:),'Type1 SHG', Lambda, Pol, refIdx, k, w);
        obj(i,j,:) = GouyPhase+squeeze(DK(i,j,:)).*Z';
        std_(i,j)  = std(squeeze(obj(i,j,:)));
        mean_(i,j) = mean(squeeze(obj(i,j,:)));
    end
end

minMatrix = min(std_(:));
[row_std,col_std] = find(std_==minMatrix);

minMatrix = min(abs(mean_(:)));
[row_mean,col_mean] = find(abs(mean_)==minMatrix);

Tpm(row_std,1,1)
Dt(1,col_std,1)/2

% [X,Y] = meshgrid(1:MatDim,1:MatDim);
% surf(X,Y,std_)