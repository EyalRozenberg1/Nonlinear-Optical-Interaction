function [ Print, figNum,  PlaneGauss_, GradType, Undepleted] = ReqPrints()

% Create temperature gradient
GradTypeVec             = {'Const','Linear', 'ApodizationMain', 'Exp_5x5x50','Exp_4x4x50', 'Exp_3x3x30', 'ArcTan', 'Apodization3','Apodization1','Apodization2','Apodization3_Old' ,'Apodization1.5','Apodization1.6','Apodization3.5','Exponential'};
GradType                = char(GradTypeVec(1));

% Define Simulation Type: Gauss=0 OR Plane wave=1
PlaneGauss_             = 0;

% Define Pump as Undepleted=1, depleted=0
Undepleted              = 0;

% set to '1' if plot is required
 Print.DeltaK           = 0;
 Print.Temperature      = 0;
 Print.RefIndex         = 0;
 Print.Amplitude        = 0; % TODO: does not support Gauss wave yet
 Print.Intensity        = 1; % TODO: THG does not support Gauss wave yet
 Print.GouyPhase        = 0;
 
 Print.P2w_vs_Temp      = 0; % NOTE: The result will be given for const gradient only (while using T.pm)
 Print.P2wVsTm          = 0;
 Print.P2w_vs_GradDiff  = 0;
 Print.P2wVsw0          = 0;
 Print.BW               = 0;
 Print.P2w_vs_Pw        = 0;
 Print.NormST           = 0; % Note: - does not support Gauss wave
                             %       - Should get positive gradient slope
 
% Gouy Phase compensation
% Print.GouyPhaseCom      = 0;
 
% Print experimental results
Print.Exp.P2wVsT        = 0;
Print.Exp.P2wVsPw       = 0;
Print.Exp.effcncyVsPw   = 0;
Print.Exp.P2w_vs_GradDiff = 0;

Print.Exp.Results       = ( Print.Exp.P2wVsPw || Print.Exp.P2wVsT || Print.Exp.effcncyVsPw || Print.Exp.P2w_vs_GradDiff);
 
% reset number of increasing figures number
figNum                  = 0;

end