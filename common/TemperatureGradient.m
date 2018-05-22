function [ TempGrad ] = TemperatureGradient( T, L, CrystalPropAxis, GradType, NumOfPoints, samples)
% the function creats a temperature gradient function %

switch GradType
    case 'Const'
        TempGrad    = zeros(1,NumOfPoints) + T.pm;
        TempGrad_tv = linspace(0, T.D_tv, samples)';
        TempGrad    = bsxfun(@minus,TempGrad, TempGrad_tv);
    case 'Linear'
        TempGrad = CrystalPropAxis.*(T.max - T.min)./L + T.min;
        TempGrad_tv = linspace(0, T.D_tv, samples)';
        TempGrad    = bsxfun(@minus,TempGrad, TempGrad_tv);
    case'ApodizationMain'
        % Three liner parts %
%         1.
%         ApodizedDeltaT  = 0.85*(T.max - T.min);
%         ApodizedLength  = 0.2*L;
%         2.
        ApodizedDeltaT  = 0.3778*(T.max - T.min);
        ApodizedLength  = 0.2*L;
        
        
        AdiabaticLength = L - ApodizedLength;
        T1              = T.min + 0.5*ApodizedDeltaT;
        T2              = T.max - 0.5*ApodizedDeltaT;
        m1              = (T1 - T.min)/(0.5*ApodizedLength);
        m2              = (T2 - T1)/AdiabaticLength;
        m3              = (T.max - T2)/(0.5*ApodizedLength);
        
        slop1 = m1*CrystalPropAxis + T.min;
        slop2 = m2*(CrystalPropAxis - 0.5*ApodizedLength) + T1;
        slop3 = m3*(CrystalPropAxis - L) + T.max;
        TempGrad = [slop1(1:round(NumOfPoints*(0.5*ApodizedLength/L)) - 1),slop2(round(NumOfPoints*(0.5*ApodizedLength/L)):round(NumOfPoints*(1 - 0.5*ApodizedLength/L)) - 1),slop3(round(NumOfPoints*(1 - 0.5*ApodizedLength/L)):end)];
%         TempGrad = smooth(TempGrad,500);
        TempGrad = TempGrad';
        
%         coefs = polyfit(CrystalPropAxis, TempGrad, 3);
%         Poly = polyval(coefs,CrystalPropAxis);
%         TempGrad = Poly;
        
case'Exp_5x5x50'
        CopperLegth     = 10e-3;
%         CopperLegth     = 4e-3;
        AdiabaticLength = L - 2*CopperLegth;

        m1              = 0;
        m2              = (T.max - T.min)/AdiabaticLength;
        m3              = 0;
        
        slop1 = m1*CrystalPropAxis + T.min;
        slop2 = m2*(CrystalPropAxis - CopperLegth) + T.min;
        slop3 = m3*(CrystalPropAxis - L) + T.max;
        
        TempGrad = [slop1(1:round(NumOfPoints*(CopperLegth/L)) - 1),slop2(round(NumOfPoints*(CopperLegth/L)):round(NumOfPoints*(1 - CopperLegth/L)) - 1),slop3(round(NumOfPoints*(1 - CopperLegth/L)):end)];
%         TempGrad = smooth(TempGrad,5777);
%         TempGrad = TempGrad';

case'Exp_4x4x50'
        CopperLegth     = 8e-3;
        AdiabaticLength = L - 2*CopperLegth;

        m1              = 0;
        m2              = (T.max - T.min)/AdiabaticLength;
        m3              = 0;
        
        slop1 = m1*CrystalPropAxis + T.min;
        slop2 = m2*(CrystalPropAxis - CopperLegth) + T.min;
        slop3 = m3*(CrystalPropAxis - L) + T.max;
        
        TempGrad = [slop1(1:round(NumOfPoints*(CopperLegth/L)) - 1),slop2(round(NumOfPoints*(CopperLegth/L)):round(NumOfPoints*(1 - CopperLegth/L)) - 1),slop3(round(NumOfPoints*(1 - CopperLegth/L)):end)];
%         TempGrad = smooth(TempGrad,5777);
%         TempGrad = TempGrad';

case'Exp_3x3x30'
        CopperLegth     = 8e-3;
        AdiabaticLength = L - 2*CopperLegth;

        m1              = 0;
        m2              = (T.max - T.min)/AdiabaticLength;
        m3              = 0;
        
        slop1 = m1*CrystalPropAxis + T.min;
        slop2 = m2*(CrystalPropAxis - CopperLegth) + T.min;
        slop3 = m3*(CrystalPropAxis - L) + T.max;
        
        TempGrad = [slop1(1:round(NumOfPoints*(CopperLegth/L)) - 1),slop2(round(NumOfPoints*(CopperLegth/L)):round(NumOfPoints*(1 - CopperLegth/L)) - 1),slop3(round(NumOfPoints*(1 - CopperLegth/L)):end)];
%         TempGrad = smooth(TempGrad,3777);
%         TempGrad = TempGrad';
    
    case 'ArcTan'
        % To reduce Gouy Phase
        xi       = 2.84;
        z0       = L/(2*xi);
        TempGrad = T.pm + 0.5 * (T.max - T.min) * atan((CrystalPropAxis-L/2)/z0)/max(atan((CrystalPropAxis-L/2)/z0));
%         TempGrad = T.pm + 1.232/1.615 * 0.5 * (T.max - T.min) * atan((CrystalPropAxis-L/2)/z0)/max(atan((CrystalPropAxis-L/2)/z0));

    case'Apodization1'
        % Three liner parts %
        slop1 = 4*CrystalPropAxis(1:round(NumOfPoints*0.1)).*(T.max - T.min)/L + T.min;
        slop2 = 0.5*CrystalPropAxis(1:round(NumOfPoints*0.8)).*(T.max - T.min)/L + slop1(end);
        slop3 = 4*CrystalPropAxis(1:round(NumOfPoints*0.1)).*(T.max - T.min)/L + slop2(end);
        TempGrad = [slop1,slop2,slop3];
    case'Apodization1.5'
        % Three liner parts %
        slop1 = 4*CrystalPropAxis(1:round(NumOfPoints*0.1)).*(T.max - T.min)/L + T.min;
        slop2 = 0.1*CrystalPropAxis(1:round(NumOfPoints*0.8)).*(T.max - T.min)/L + slop1(end);
        slop3 = 4*CrystalPropAxis(1:round(NumOfPoints*0.1)).*(T.max - T.min)/L + slop2(end);
        TempGrad = [slop1,slop2,slop3];
    case'Apodization1.6'
        % Three liner parts %
        slop1 = 2*CrystalPropAxis(1:round(NumOfPoints*0.3)).*(T.max - T.min)/L + T.min;
        slop2 = 0.1*CrystalPropAxis(1:round(NumOfPoints*0.4)).*(T.max - T.min)/L + slop1(end);
        slop3 = 2*CrystalPropAxis(1:round(NumOfPoints*0.3)).*(T.max - T.min)/L + slop2(end);
        TempGrad = [slop1,slop2,slop3];
    case'Apodization2'
        % A + Bx^3 %
        A = (T.max+T.min)/2;
        B = (T.max-T.min)/(2*(L/2)^3);
        TempGrad = A + B*(CrystalPropAxis-L/2).^3;
    case'Apodization3_Old'
        % NOTE: This one incrase from the Tmax value 
        C = 400;
        A = (T.max+T.min)/2;
        B = (T.max-T.min)/2*(L/2)^-9 - C*(L/2);
        TempGrad = A + B*(CrystalPropAxis-L/2).^9 + C*(CrystalPropAxis);
    case'Apodization3'
        % A + Bx^9 + Cx %
        C = 500; % C is will is in cherge of the relation of the slope in the middle and the edges
        A = (T.max+T.min)/2;
        B = ((T.max-T.min)/2 - C*(L/2))*(L/2)^-9;
        TempGrad = A + B*(CrystalPropAxis-L/2).^9 + C*(CrystalPropAxis-L/2);
        plot(TempGrad)
     case'Apodization3.5'
        % A + Bx^9 + Cx %
        C = 20; % C is will is in cherge of the relation of the slope in the middle and the edges
        A = (T.max+T.min)/2;
        B = ((T.max-T.min)/2 - C*(L/2))*(L/2)^-9;
        TempGrad = A + B*(CrystalPropAxis-L/2).^9 + C*(CrystalPropAxis-L/2);
        plot(TempGrad)
    case 'Exponential'
        Arg         = log(T.max)/L;
        TempGrad    = exp(Arg*CrystalPropAxis) + T.min;
end

