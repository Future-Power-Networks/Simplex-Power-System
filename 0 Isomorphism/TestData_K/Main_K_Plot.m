% This function analyzes the xi

clear all
clc
close all

mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

%% Enables
Enable_SaveFigure           = 1;
Enable_NonZeroXi            = 0;
Enable_ClearNodeNumber      = 0;
Enable_ParticipationColor   = 0;
Enable_NodeImpedanceColor   = 0;

FiedlerAbsMax = 0.7;

%% Load data
% DataName = 'K_68Bus_IBR_Load_Data';
% DataName = 'K_68Bus_IBR_Data';
% DataName = 'K_68Bus_IBR_17_Data';
% DataName = 'K_68Bus_IBR_17_14_Data';
DataName = 'K_68Bus_IBR_17_14_7_Data';

Data = load(DataName).SaveData;

%%
KH              = Data.KH;
YbusVI          = Data.YbusVI;
YbusVIF         = Data.YbusVIF;
GbusVI          = Data.GbusVI;
GbusVIF         = Data.GbusVIF;
YbusOrigin      = Data.YbusOrigin;
Index_Vbus      = Data.Index_Vbus;
Index_Ibus      = Data.Index_Ibus;
Index_Fbus      = Data.Index_Fbus;
Index_Ebus      = Data.Index_Ebus;
Order_Old2New   = Data.Order_Old2New;
Order_New2Old   = Data.Order_New2Old;

%%
FigNum = 0;
ColorRGB();
FigSize = [0.1 0.1 0.5 0.75];

FigNum = FigNum + 1;
figure(FigNum)
set(gcf,'units','normalized','outerposition',FigSize);

%% Plot graph
GraphMatrix = NormMatrixElement(YbusOrigin,'DiagFlag',0);
GraphData = graph(GraphMatrix,'upper');
GraphFigure = plot(GraphData); grid on; hold on;
highlight(GraphFigure,GraphData,'EdgeColor',[0,0,0],'LineWidth',1.1);       % Change all edges to black by default
highlight(GraphFigure,GraphData,'NodeColor',[0,0,0]);                    	% Change all nodes to black by default
highlight(GraphFigure,GraphData,'MarkerSize',5);

%% Set voltage node
highlight(GraphFigure,Index_Vbus,'NodeColor',[0,0,0]);      	% Change all voltage node to black by default

%% Set floating node
highlight(GraphFigure,Index_Fbus,'NodeColor',[0.7,0.7,0.7]);   	% Change all floating node to gray by default

%% Reduce pure empty node
highlight(GraphFigure,Index_Ebus,'MarkerSize',1);  	% Reduce the size of empty node, i.e., node without apparatus and passive load
                                                    % Empty node is floating node for system without any passive load
% highlight(GraphFigure,Index_Fbus,'Marker','o');

%% Calculation
[Phi,Xi,Psi] = eig(KH);
PhiInv = inv(Phi);
Xi = diag(Xi);
[~,Index_XiMin] = min(real(Xi));
XiMin = Xi(Index_XiMin);
if Enable_NonZeroXi
    if abs(XiMin)<=1e-4                    % Check if xi_min is zero
        Xi_ = Xi;
        Xi_(Index_XiMin) = inf;
        [~,Index_XiMin] = min(real(Xi_));
        XiMin = Xi(Index_XiMin);
    end
end

PhiRightMin = Phi(:,Index_XiMin);
PhiLeftMin = transpose(PhiInv(Index_XiMin,:));
FiedlerVec = PhiRightMin.*PhiLeftMin;
FiedlerAbs = abs(FiedlerVec);
% FiedlerAbsMax = max(FiedlerAbs);
FiedlerAbsMax = max([max(FiedlerAbs),FiedlerAbsMax]);
FiedlerAbsNorm = FiedlerAbs/FiedlerAbsMax;

%% Color bar
ColorStepSize = 100;
% ColorLower = [1,0.95,0];
% ColorUpper = [1,0,0];
ColorLower = [0,1,1];       % Light blue
ColorUpper = [0,0,1];       % Dark blue
GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStepSize)';
GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStepSize)';
GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStepSize)';
% colormap([GradRed GradGreen GradBlue]);
% colorbar();
% caxis([0 FiedlerAbsMax]);

%% Set current node
highlight(GraphFigure,Index_Ibus,'NodeColor',[0,1,0]);  % Change all current node to green by default

if Enable_ParticipationColor
ColorFactor = FiedlerAbsNorm((length(Index_Vbus)+1):end);
for k = 1:length(Index_Ibus)
    NodeColor{k} = (ColorUpper - ColorLower) * ColorFactor(k) + ColorLower;
    highlight(GraphFigure,Index_Ibus(k),'NodeColor',NodeColor{k});
end
end

%% Clear all node number
if Enable_ClearNodeNumber
    for k = 1:length(Order_Old2New)
        GraphFigure.NodeLabel{k} = '';
    end
end

%% Calculate impedance
GbusVIF_Diag = diag(GbusVIF);
Znode = GbusVIF_Diag(Order_New2Old);     % Convert the order back to its origin
% GraphFigure.NodeFontSize = 7;

% Convert admittance to impedance at voltage node
for k = 1:length(Index_Vbus)
    Znode(Index_Vbus(k)) = 1/Znode(Index_Vbus(k));        
end

% Add color for current node based on impedance's phase
if Enable_NodeImpedanceColor
    for k = 1:length(Index_Ibus) 
        Zimag = imag(Znode(Index_Ibus(k)));
        if Zimag>=0
            highlight(GraphFigure,Index_Ibus(k),'NodeColor',[1,1,0]);   % Yellow for inductive
        else
            highlight(GraphFigure,Index_Ibus(k),'NodeColor',[0,1,0]);   % Green for capacitive
        end
    end
end

ZnodeAbs = abs(Znode);

XData = GraphFigure.XData';
YData = GraphFigure.YData';
ZData = ZnodeAbs;

XData = XData(Order_Old2New);
YData = YData(Order_Old2New);
ZData = ZData(Order_Old2New);

% Remove the voltage bus impedance to smooth the map
XData = XData(length(Index_Vbus)+1:end,:);
YData = YData(length(Index_Vbus)+1:end,:);
ZData = ZData(length(Index_Vbus)+1:end,:);

%% Plot impedance value to node
if 0
    for k = 1:length(Order_Old2New)
        GraphFigure.NodeLabel{k} = num2str(ZnodeAbs(k),2);
    end
end

%% Plot impedance heat map
% Set the value at four cornors and four edges
% XData_1 = linspace(-4.5,4.5,10)';
% YData_1 = ones(size(XData_1))*(4.5);
% ZData_1 = zeros(size(XData_1));
% XData_2 = XData_1;
% YData_2 = -YData_1;
% ZData_2 = ZData_1;
% XData_3 = YData_1;
% YData_3 = XData_1;
% ZData_3 = ZData_1;
% XData_4 = -YData_1;
% YData_4 = XData_1;
% ZData_4 = ZData_1;
% XData = [XData; XData_1; XData_2; XData_3; XData_4];
% YData = [YData; YData_1; YData_2; YData_3; YData_4];
% ZData = [ZData; ZData_1; ZData_2; ZData_3; ZData_4];

% Plot heat map
PlotHeatMap(XData,YData,ZData,1,[0,0.25]);

% Get the max
ZDataMax = max(ZData)

% Move graph to top
uistack(GraphFigure,'top');

%% Set Figure Lim
FigureMargin = 0.3;
% xmax = max(abs(XData));
% ymax = max(abs(YData));
% xymax = max(xmax,ymax);
% xlim([-xymax-FigureMargin,xymax+FigureMargin]);
% ylim([-xymax-FigureMargin,xymax+FigureMargin]);
xlim([min(GraphFigure.XData)-FigureMargin,max(GraphFigure.XData)+FigureMargin]);
ylim([min(GraphFigure.YData)-FigureMargin,max(GraphFigure.YData)+FigureMargin]);

%% Save
if Enable_SaveFigure
    print(figure(1),['Graph_' DataName '.png'],'-dpng','-r600');
end