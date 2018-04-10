function FormatPlotFontSizeNameLine(h_ax,h_fig,font_size,font_name,line_width)

title_factor=1; %4/3;
set(h_ax, 'FontName', font_name);
set(h_ax,'Linewidth',2);
set(h_ax,'Fontsize',font_size,'Fontweight','bold');
set(h_ax,'Box','on');
set(h_fig,'Color',[1,1,1]);
set(get(h_ax,'XLabel'),'Fontsize',ceil(font_size*title_factor),'Fontweight','bold');
set(get(h_ax,'YLabel'),'Fontsize',ceil(font_size*title_factor),'Fontweight','bold');
set(get(h_ax,'Title'),'Fontsize',ceil(font_size*title_factor),'Fontweight','bold');
% set(get(h_ax,'Legend'),'Fontsize',ceil(font_size/2),'Fontweight','bold');
hlist=get(h_ax,'Children');
for index=1:length(hlist)
    set(hlist(index),'LineWidth',line_width);
end

%anat
% set(findall(gcf,'type','text'),'FontSize',font_size)
% set(gca,'FontSize',font_size)

set(findall(gcf,'type','text'),'FontName',font_name)
set(gca,'FontName',font_name)

