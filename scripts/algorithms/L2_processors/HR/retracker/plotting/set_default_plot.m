% close all;

global defaultLineLineWidth defaultLineMarkerSize

defaultLineLineWidth=3;
defaultLineMarkerSize=5;

textbox_fontsize=24;
colorbar_fontsize=24;
legend_fontsize=24;

% width = 8;
% height = 6;
VISIBLE=0;

% load('.\plotting\colormapDEM.mat');
% colormapATDD = colormap('hot');
% colormapATDD(1,:)=[1 1 1];
set(0,'defaultLineLineWidth',defaultLineLineWidth);   % set the default line width
set(0,'defaultLineMarkerSize',defaultLineMarkerSize);  % set the default line marker size
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',24);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',24);

% Make ALL the plots invisible
if VISIBLE == 0
    set(0,'defaultFigureVisible','off');
else
    set(0,'defaultFigureVisible','on');
end

% Set the default Size for display
mida = get(0,'ScreenSize');
mida(3:4)=[1536, 864]; % [1920,1080];
% set(0,'defaultFigurePosition', [defpos(1) defpos(2)-50 width*100, height*100]);
set(0,'defaultFigurePosition',mida);
close

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','centimeters');
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
set(0,'defaultFigurePaperUnits','points');
set(0,'defaultFigurePaperPosition', mida);

% plot(...)
% print('figName','-r100','-dpng')