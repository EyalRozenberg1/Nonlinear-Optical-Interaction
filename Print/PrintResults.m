function [figNum] = PrintResults(PinPeak, PlaneGauss_, figNum, CrystalPropAxis, DeltaK, Efficiency, Lambda, T, A, P, n, I, TempGrad, InteractionType, Print, NumOfPoints, printBW, I_from_A, BeamWidthVec, xi)

% this function will print required results
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultAxesFontName','times')
set(0,'DefaultTextFontSize',18)
set(0,'DefaultTextFontName','times')

if (printBW)
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    
    plot(Lambda.inVec, Efficiency,'b*-','Linewidth',2);
    xlim([min(Lambda.inVec) max(Lambda.inVec)]);
    
    xlabel('\lambda_p_u_m_p [m]');
    ylabel('Conversion Efficiency');
    
    title(['Conversion Efficiency vs \lambda for ',InteractionType,' P_I_n= ',num2str(PinPeak),' W']);
%     [M,S] = max(Efficiency);
    [~,S] = min(abs(Lambda.inVec-Lambda.in1));
    text(Lambda.inVec(S),Efficiency(S), {['\leftarrow \eta= ', num2str(Efficiency(S)),' ,\lambda= ', num2str(Lambda.inVec(S))]});

elseif(Print.P2w_vs_Temp)
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    % set(FigHandle, 'Position', [83 374 576 512]);

    plot(T,Efficiency,'b*-','Linewidth',2);
    ylabel('Conversion Efficiency');
    str = sprintf('Conversion Efficiency Vs T for %s, P_I_n= %s W',InteractionType, num2str(PinPeak));
    title(str);
    xlabel('T [^\circC]');
    xlim([T(1) T(end)])
    [M,S] = max(Efficiency);
    text(T(S), Efficiency(S), {['\leftarrow \eta= ', num2str(M),' ,T = ',num2str(T(S))] ['  \DeltaK\bulletL = ',num2str(DeltaK(S))]});
    
    % set data cursor
    dcm_obj = datacursormode(FigHandle);
    set(dcm_obj,'UpdateFcn',{
        @myupdatefcn2, 	...
        DeltaK});
    
    
elseif(Print.P2w_vs_GradDiff)
	figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    
    plot(2*T, Efficiency,'b*-','Linewidth',2);
    xlim([min(2*T) max(2*T)]);
    
    xlabel('\DeltaT [^\circC]');
    ylabel('Conversion Efficiency');

    title(['Conversion Efficiency Vs \DeltaT for ',InteractionType,', P_I_n= ',num2str(PinPeak),' W']);

    [M,S] = max(Efficiency);
    text(2*T(S), Efficiency(S), {['\leftarrow \eta= ', num2str(M),' ,\DeltaT= ', num2str(2*T(S))]});


elseif(Print.P2wVsTm)
	figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    % set(FigHandle, 'Position', [83 374 576 512]);

    plot(T, Efficiency,'b*-','Linewidth',2);
    xlim([min(T) max(T)]);
    
    xlabel('Temp middle [^\circC]');
    ylabel('Conversion Efficiency');

    title(['Conversion Efficiency Vs Tm for ',InteractionType,', P_I_n= ',num2str(PinPeak),' W']);

    [M,S] = max(Efficiency);
    text(T(S), Efficiency(S), {['\leftarrow \eta= ', num2str(M),' ,Tm= ',num2str(T(S))]});

    
elseif(Print.P2wVsw0)
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    % set(FigHandle, 'Position', [83 374 576 512]);
    plot(BeamWidthVec,Efficiency,'b*-','Linewidth',2);
%     str = sprintf('Conversion Efficiency Vs Beam waist width for %s, P_I_n= %s W',InteractionType, num2str(PinPeak));
    title(['Conversion Efficiency Vs Beam waist width for ',InteractionType,' P_I_n= ', num2str(PinPeak),' W']);
    ylabel('Conversion Efficiency');
    xlabel('w_0 [m]');
    [M,S] = max(Efficiency);
    text(BeamWidthVec(S), Efficiency(S), ['\leftarrow \eta= ', num2str(M),' ,w_o= ',num2str(BeamWidthVec(S))]);
    
    
elseif(Print.P2w_vs_Pw)
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    % set(FigHandle, 'Position', [83 374 576 512]);
    % print P2w VS Pw
    if(PlaneGauss_)
        P1  = P.in1/10^4; % convertion from W/m^2 to W/cm^2
        P2  = P.out/10^4; % convertion from W/m^2 to W/cm^2
        plot(P1,P2,'k*-','Linewidth',2);
        str = sprintf('I2w Vs Iw for %s',InteractionType);
        title(str);
        xlabel('I\omega [W/cm^2]');
        ylabel('I2\omega [W/cm^2]');
        hold on;
        coefs = polyfit(P1, P2, 2);
        P2wVsPwPoly = polyval(coefs,0:1:max(P1));
        plot(0:1:max(P1),P2wVsPwPoly,'m','Linewidth',0.5);
        a1 = coefs(2);
        a2 = coefs(1);
        polyfit_str = ['P2\omega = ' num2str(a1) ' P\omega + ' num2str(a2) ' P\omega^2'];
        %text(7,7,polyfit_str)
        hold off;
        legend('Simulation',polyfit_str);
    else
        plot(P.in1,P.out,'k*-','Linewidth',2);
        str = sprintf('P2w Vs Pw for %s',InteractionType);
        title(str);
        xlabel('P\omega [W]');
        ylabel('P2\omega [W]');
        hold on;
        coefs = polyfit(P.in1, P.out, 2);
        P2wVsPwPoly = polyval(coefs,0:1:max(P.in1));
        plot(0:1:max(P.in1),P2wVsPwPoly,'m','Linewidth',0.5);
        a1 = coefs(2);
        a2 = coefs(1);
        polyfit_str = ['P2\omega = ' num2str(a1) ' P\omega + ' num2str(a2) ' P\omega^2'];
        %text(7,7,polyfit_str)
        hold off;
        legend('Simulation',polyfit_str);
    end
    figNum = figNum + 1;
    FigHandle = figure(figNum);set(gcf,'color','white');
    % set(FigHandle, 'Position', [83 374 576 512]);
    
    % print Conversion Efficiency - P2w/Pw Vs Pw
    ConvEfficiency = 100*P.out./P.in1;
    ConvEfficiency(find(isnan(ConvEfficiency) == 1))=0;
    if(PlaneGauss_)
        P1  = P.in1/10^4; % convertion from W/m^2 to W/cm^2 
        plot(P1,ConvEfficiency,'k*-','Linewidth',2);
        str = sprintf('Conversion Efficiency Vs Iw for %s',InteractionType);
        title(str);
        xlabel('I\omega [W/cm^2]');
        ylabel('% I2\omega/I\omega');
        hold on;
        coefs2 = polyfit(P1, ConvEfficiency, 1);
        PrcntgPoly = polyval(coefs2,0:1:max(P1));
        plot(0:1:max(P1),PrcntgPoly,'m','Linewidth',0.5);
        b1 = coefs2(1);
        polyfit_str2 = ['P2\omega/P\omega = ' num2str(b1) ' P\omega'];
        %text(7,7,polyfit_str2)
        hold off;
        legend('Simulation',polyfit_str2);
    else
        plot(P.in1,ConvEfficiency,'k*-','Linewidth',2);
        str = sprintf('Conversion Efficiency Vs Pw for %s',InteractionType);
        title(str);
        xlabel('P\omega [W]');
        ylabel('% P2\omega/P\omega');
        hold on;
        coefs2 = polyfit(P.in1, ConvEfficiency, 1);
        PrcntgPoly = polyval(coefs2,0:1:max(P.in1));
        plot(0:1:max(P.in1),PrcntgPoly,'m','Linewidth',0.5);
        b1 = coefs2(1);
        polyfit_str2 = ['P2\omega/P\omega = ' num2str(b1) ' P\omega'];
        %text(7,7,polyfit_str2)
        hold off;
        legend('Simulation',polyfit_str2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  short simulation graphs %%%%%%%%%%%%%%%%%%%%
else
    if (Print.DeltaK)
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % print Delta K as a function of Crystal length
        plot(CrystalPropAxis, DeltaK,'LineWidth',2);
        title('\DeltaK as a function of crystal length');
        xlabel('Crystal Interaction Length [m]');
        ylabel('\DeltaK [1/m]');
        if(DeltaK(1) ~= DeltaK(end))
            text(CrystalPropAxis(1), DeltaK(1), ['\leftarrow T=', num2str(T.min)]);
            text(CrystalPropAxis(end), DeltaK(end), ['\leftarrow T=', num2str(T.max)]);
            if(T.max>=T.min)
                text(CrystalPropAxis(min(find(DeltaK>=0))), DeltaK(min(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(min(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(min(find(DeltaK>=0)))),' [m]']});
            else
                text(CrystalPropAxis(max(find(DeltaK>=0))), DeltaK(max(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(max(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(max(find(DeltaK>=0)))),' [m]']});
            end
        else
            text(CrystalPropAxis(floor(NumOfPoints/1.5)), DeltaK(floor(NumOfPoints/1.5)), ['T=', num2str(T.pm)]);
            text(CrystalPropAxis(floor(NumOfPoints/4)), DeltaK(floor(NumOfPoints/4)), ['\DeltaK=', num2str(DeltaK(floor(NumOfPoints/2)))]);
        end
    end
    
    
    if (Print.Temperature)
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % print Temperature gradient as a function of Crystal length
        plot(CrystalPropAxis, TempGrad);
        title('Temperature gradient as a function of crystal length');
        xlabel('Crystal Interaction Length [m]');
        ylabel('T [^\circC]');
        if(DeltaK(1) ~= DeltaK(end))
            if(T.max>=T.min)
                text(CrystalPropAxis(min(find(DeltaK>=0))), TempGrad(min(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(min(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(min(find(DeltaK>=0)))),' [m]']});
            else
                text(CrystalPropAxis(max(find(DeltaK>=0))), TempGrad(max(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(max(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(max(find(DeltaK>=0)))),' [m]']});
            end
        else
            text(CrystalPropAxis(floor(NumOfPoints/2)), TempGrad(1), ['T=', num2str(T.pm)]);
        end
    end
        
    if (Print.GouyPhase)
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % calculate Gouy Phase
        L         = CrystalPropAxis(end);
        z0        = L/(2*xi);
        GouyPhase = -atan((CrystalPropAxis-L/2)/z0).';
        plot(CrystalPropAxis, GouyPhase, 'r', 'Linewidth', 2);
        title(['Gouy Phase for \xi = ',num2str(xi),', \DeltaK\bulletL = ',num2str(L*(DeltaK(end)-DeltaK(1)))]);
        xlabel('Crystal Interaction Length [m]');
        ylabel('Gouy Phase');
        hold on;
        plot(CrystalPropAxis,DeltaK.*CrystalPropAxis','b', 'Linewidth', 2);
        plot(CrystalPropAxis,GouyPhase+DeltaK.*CrystalPropAxis','k', 'Linewidth', 2);
        legend('Gouy Phase', '\DeltaK\bulletL','Gouy Phase+\DeltaK\bulletL');
        hold off;
        if(DeltaK(1) ~= DeltaK(end))
            text(CrystalPropAxis(1), DeltaK(1)*L, ['\leftarrow T=', num2str(T.min)]);
            text(CrystalPropAxis(end), DeltaK(end)*L, ['\leftarrow T=', num2str(T.max)]);
            if(T.max>=T.min)
                text(CrystalPropAxis(min(find(DeltaK>=0))), DeltaK(min(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(min(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(min(find(DeltaK>=0)))),' [m]']});
            else
                text(CrystalPropAxis(max(find(DeltaK>=0))), DeltaK(max(find(DeltaK>=0))), {['   \leftarrow PM: T=',num2str(TempGrad(max(find(DeltaK>=0)))),' [^\circC]'],['               Pos=',num2str(CrystalPropAxis(max(find(DeltaK>=0)))),' [m]']});
            end
        end
    end
    
    if (Print.RefIndex)
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % print Refractive index as a function of Crystal length
        if strcmp('Type1 SHG',InteractionType)
            plot(CrystalPropAxis, n.in1,'r', 'Linewidth', 2); hold on;
            plot(CrystalPropAxis, n.out,'b', 'Linewidth', 2); hold off;
            title({'Refractive index as a function of Crystal Crystal length','Type 1 oo(z)->e(y)'});
            legend('no(z)','ne(y)');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Refractive Index');
            
            figNum = figNum + 1;
            FigHandle = figure(figNum);set(gcf,'color','white');
            plot(CrystalPropAxis, n.in1 - n.out,'k', 'Linewidth', 2);
            title({'Refractive index difference as a function of Crystal Crystal length','Type 1 oo(z)->e(y)'});
            legend('\Deltan');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Refractive Index');
        end
        
        if strcmp('Type1 THG',InteractionType)
            plot(CrystalPropAxis, 0.5*(n.in1+n.in2), CrystalPropAxis, n.out);
            title({'Refractive index as a function of Crystal Crystal length','Type 1 oo(z)->e(y)'});
            legend('(no_1(z)+no_2(z))/2','ne(y)');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Refractive Index');
        end
    end
    
    
    if(Print.Amplitude) % TODO: does not support Gauss wave yet
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % print Amplitude as a function of Crystal length
        if strcmp('Type1 SHG',InteractionType)
            AinPower2  = A.in1.*conj(A.in1);
            AoutPower2 = A.out.*conj(A.out);
            plot(CrystalPropAxis,[AinPower2(:,1) AoutPower2(:,1)]);
            title('Amplitude as a function of crystal length');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Amplitude [Ws/m^2F]');
            legend('|Ain|^2(\omega)','|Aout|^2(2\omega)');
        end
        if strcmp('Type1 THG',InteractionType)
            Ain1Power2  = A.in1.*conj(A.in1);
            Ain2Power2  = A.in2.*conj(A.in2);
            AoutPower2 = A.out.*conj(A.out);
            plot(CrystalPropAxis,[Ain1Power2(:,1) Ain2Power2(:,1) AoutPower2(:,1)]);
            title('Amplitude as a function of crystal length');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Amplitude [Ws/m^2F]');
            legend('|Ain_1|^2(\omega)','|Ain_2|^2(\omega)','|Aout|^2(\omega)');
        end
    end
    
    
    if(Print.Intensity)
        figNum = figNum + 1;
        FigHandle = figure(figNum);set(gcf,'color','white');
        % set(FigHandle, 'Position', [83 374 576 512]);
        % print intensity as a function of Crystal length
        if strcmp('Type1 SHG',InteractionType)
            %figure;
            if (PlaneGauss_)
                AinPower2  = A.in1.*conj(A.in1); AinPower2  = AinPower2.';
                AoutPower2 = A.out.*conj(A.out); AoutPower2 = AoutPower2.';
                IntensityIn  = I_from_A(mean(n.in1),AinPower2)/10^4; % convertion from W/m^2 to W/cm^2;
                IntensityOut = I_from_A(mean(n.out),AoutPower2)/10^4; % convertion from W/m^2 to W/cm^2;
                plot(CrystalPropAxis,IntensityIn,'r','LineWidth',2);hold on;
                plot(CrystalPropAxis,IntensityOut,'b','LineWidth',2);
                plot(CrystalPropAxis,IntensityIn+IntensityOut,'k','LineWidth',2); hold off;
                str = sprintf('Intensity as a function of crystal length, I_I_n= %s [GW/cm^2]', num2str((I.in1/10^4)/1e9));
                title(str);
                %title('Intensity as a function of crystal length');
                xlabel('Crystal Interaction Length [m]');
                ylabel('Intensity [W/cm^2]');
                legend('Iin(\omega)','Iout(2\omega)','Iin(\omega)+Iout(2\omega)');
                if(T.max>=T.min)
                    text(CrystalPropAxis(min(find(DeltaK>=0))), IntensityOut(min(find(DeltaK>=0))), '\leftarrow PM');
                else
                    text(CrystalPropAxis(max(find(DeltaK>=0))), IntensityOut(max(find(DeltaK>=0))), '\leftarrow PM');
                end
            else
                plot(CrystalPropAxis,P.in1,'r','LineWidth',2); hold on;
                plot(CrystalPropAxis,P.out,'b','LineWidth',2);
                plot(CrystalPropAxis,P.out+P.in1,'k','LineWidth',2); hold off;
                str = sprintf('Power as a function of crystal length, P_I_n= %s W', num2str(PinPeak));
                title(str);
                %title('Power as a function of crystal length');
                xlabel('Crystal Interaction Length [m]');
                ylabel('Peak Power [W]');
                legend('Pin(\omega)','Pout(2\omega)', 'Pin(\omega)+Pout(2\omega)');
                if(T.max>=T.min)
                    text(CrystalPropAxis(min(find(DeltaK>=0))), P.out(min(find(DeltaK>=0))), '\leftarrow PM');
                else
                   text(CrystalPropAxis(max(find(DeltaK>=0))), P.out(max(find(DeltaK>=0))), '\leftarrow PM'); 
                end
            end
        end % TODO: does not support Gauss wave yet
        if strcmp('Type1 THG',InteractionType)
            Ain1Power2  = A.in1.*conj(A.in1); Ain1Power2  = Ain1Power2(:,1).';
            Ain2Power2  = A.in2.*conj(A.in2); Ain2Power2  = Ain2Power2(:,1).';
            AoutPower2  = A.out.*conj(A.out); AoutPower2  = AoutPower2(:,1).';
            IntensityIn1  = I_from_A(n.in1,Ain1Power2)/10^4; % convertion from W/m^2 to W/cm^2;
            IntensityIn2  = I_from_A(n.in2,Ain2Power2)/10^4; % convertion from W/m^2 to W/cm^2;
            IntensityOut  = I_from_A(n.out,AoutPower2)/10^4; % convertion from W/m^2 to W/cm^2;
            plot(CrystalPropAxis,IntensityIn1,CrystalPropAxis,IntensityIn2,CrystalPropAxis,IntensityOut);
            title('Intensity as a function of crystal length');
            xlabel('Crystal Interaction Length [m]');
            ylabel('Intensity [W/cm^2]');
            legend('Iin_1(\omega)','Iin_2(\omega)','Iout(\omega)');
            str1 = num2str(Lambda.in1);str2 = num2str(Lambda.in2);str3 = num2str(Lambda.out);
            text(CrystalPropAxis(floor(NumOfPoints)), IntensityIn1(floor(NumOfPoints)), ['\lambda=',str1]);
            text(CrystalPropAxis(floor(NumOfPoints)), IntensityIn2(floor(NumOfPoints)), ['\lambda=',str2]);
            text(CrystalPropAxis(floor(NumOfPoints)), IntensityOut(floor(NumOfPoints)), ['\lambda=',str3]);
            if(T.max>=T.min)
                text(CrystalPropAxis(min(find(DeltaK>=0))), IntensityOut(min(find(DeltaK>=0))), '\leftarrow PM');
            else
                text(CrystalPropAxis(max(find(DeltaK>=0))), IntensityOut(max(find(DeltaK>=0))), '\leftarrow PM');
            end
        end
    end
%     if(Print.effcncy_vs_Axs)
%         figNum = figNum + 1;
%         FigHandle = figure(figNum);set(gcf,'color','white');
%         % set(FigHandle, 'Position', [83 374 576 512]);
%         % print Efficiency VS Crystal Interaction Length
%         if (PlaneGauss_)
%             AinPower2  = A.in1.*conj(A.in1); AinPower2  = AinPower2(:,1).';
%             AoutPower2 = A.out.*conj(A.out); AoutPower2 = AoutPower2(:,1).';
%             P2wVsPw = AoutPower2./AinPower2;
%         else
%             P2wVsPw = P.out./P.in1;
%         end
%         plot(CrystalPropAxis,P2wVsPw);
%         title('Efficiency Vs Crystal Interaction Length');
%         xlabel('Crystal Interaction Length [m]');
%         ylabel('$$\frac{P2\omega}{P\omega}$$','Interpreter','latex')
%         text(CrystalPropAxis(min(find(DeltaK>=0))), P2wVsPw(min(find(DeltaK>=0))), '\leftarrow PM');
%     end
end


    function txt = myupdatefcn2(~,event_obj,DeltaK_L)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        loc   = get(event_obj,'DataIndex');
        txt = {
            ['Efficieny: ',num2str(pos(2))]                     ,...
            ['T: ',num2str(pos(1))]                     ,...
            ['DeltaK*L: ',num2str(DeltaK_L(loc))]};
    end % function myupdatefcn2

function txt = myupdatefcn3(~,event_obj,DeltaK_L)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        loc   = get(event_obj,'DataIndex');
        txt = {
            ['Efficieny: ',num2str(pos(2))]                     ,...
            ['DeltaT: ',num2str(pos(1))]             ,...
            ['DeltaK*L: ',num2str(DeltaK_L(loc))]};
    end % function myupdatefcn3

end