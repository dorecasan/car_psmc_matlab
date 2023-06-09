function createfigure1(XData1, YData1, YData2, YData3)
%CREATEFIGURE1(XData1, YData1, YData2, YData3)
%  XDATA1:  line xdata
%  YDATA1:  line ydata
%  YDATA2:  line ydata
%  YDATA3:  line ydata

%  Auto-generated by MATLAB on 22-Jul-2022 16:41:10

% Create figure
figure1 = figure('Tag','ScopePrintToFigure','Color',[1 1 1],...
    'OuterPosition',[382.6 168.2 687.2 473.6]);

% Create axes
axes1 = axes('Tag','DisplayAxes1_RealMag','Parent',figure1,...
    'ColorOrder',[0 0.447058823529412 0.741176470588235;0.850980392156863 0.325490196078431 0.0980392156862745;0.929411764705882 0.694117647058824 0.125490196078431;0.494117647058824 0.184313725490196 0.556862745098039;0.466666666666667 0.674509803921569 0.188235294117647;0.301960784313725 0.745098039215686 0.933333333333333;0.635294117647059 0.0784313725490196 0.184313725490196],...
    'Position',[0.0737217598097503 0.15625 0.83127824019025 0.768749999999999]);
hold(axes1,'on');

% Create hgtransform
hgtransform('HitTest','off','Parent',axes1,...
    'Matrix',[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]);

% Create hgtransform
hgtransform('HitTest','off','Parent',axes1,...
    'Matrix',[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]);

% Create hgtransform
hgtransform('HitTest','off','Parent',axes1,...
    'Matrix',[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]);
line(XData1,YData1,'DisplayName','Demux3/4','Tag','DisplayLine2',...
    'Parent',axes1,...
    'LineWidth',3,...
    'Color',[1 0 0]);

% Create line
line(XData1,YData2,'DisplayName','Demux3/3','Tag','DisplayLine1',...
    'Parent',axes1,...
    'LineWidth',3,...
    'LineStyle','--',...
    'Color',[0 0.447058823529412 0.741176470588235]);

% Create line

% Create line
line(XData1,YData3,'DisplayName','Gain8','Tag','DisplayLine3',...
    'Parent',axes1,...
    'LineWidth',3,...
    'LineStyle','-.',...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647]);


% Create ylabel
ylabel('$\dot{e}_{y}$','FontName','Times New Roman','Interpreter','latex');

% Create xlabel
xlabel('Time (s) ','FontName','Times New Roman','Interpreter','latex');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.145454545454546 20.1454545454545]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.969281179544658 4.03071882045534]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'ClippingStyle','rectangle','FontName','Times New Roman',...
    'FontSize',14,'GridAlpha',0.4,'GridColor',[0 0 0],'GridLineStyle','--',...
    'TickLabelInterpreter','none','XGrid','on','YGrid','on');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Units','pixels',...
    'Position',[500.98438589172 303.333332274755 107.800001354218 53.0000011444092],...
    'Interpreter','none',...
    'FontSize',10,...
    'EdgeColor',[0 0 0]);

