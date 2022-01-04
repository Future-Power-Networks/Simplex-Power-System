% Author(s): Yitong Li, Yunjie Gu

clear all; 
clc;
close all;

fprintf('==================================\n')
fprintf('Start.\n')
fprintf('==================================\n')

% ### Set the name of netlist
% Other possible function: readmatrix, csvread ...
% Name_Netlist = 'netlist_Paper_SingleSGInfiniteBus.xlsx';
% Name_Netlist = 'netlist_Paper_SingleVSIInfiniteBus.xlsx';
% Name_Netlist = 'netlist_Paper_IEEE14Bus_v7.xlsx';
Name_Netlist = 'netlist_Paper_IEEE14Bus_InfBus.xlsx';

%%
% ==================================================
% Load customized data
% ==================================================
fprintf('Load customized data.\n')

% #### Load data
ListBus = xlsread(Name_Netlist,1);     
ListDevice = xlsread(Name_Netlist,2);
ListBasic = xlsread(Name_Netlist,3);
ListLine = xlsread(Name_Netlist,4);
ListLineIEEE = xlsread(Name_Netlist,5);
ListAdvance = xlsread(Name_Netlist,6);

% ### Re-arrange advanced settings
Flag_PowerFlow                  = ListAdvance(5);
Enable_CreateSimulinkModel      = ListAdvance(6);
Enable_PlotPole                 = ListAdvance(7);
Enable_PlotAdmittance           = ListAdvance(8);
Enable_PrintOutput              = ListAdvance(9);
Enable_PlotSwing                = ListAdvance(10);

% ### Re-arrange the simulation data
Fs = ListBasic(1);
Ts = 1/Fs;              % (s), sampling period
Fbase = ListBasic(2);   % (Hz), base frequency
Sbase = ListBasic(3);   % (VA), base power
Vbase = ListBasic(4);   % (V), base voltage
Ibase = Sbase/Vbase;
Zbase = Sbase/Ibase;
Ybase = 1/Zbase;
Wbase = Fbase*2*pi;

% ### Re-arrange the netlist and check error
[ListLine] = RearrangeNetlist_IEEE2Toolbox(ListLine,ListLineIEEE);
[ListBus,ListLine,ListDevice,N_Bus,N_Branch,N_Device] = RearrangeNetlist(ListBus,ListLine,ListDevice);

% ### Re-arrange the device data
[DeviceType,Para] = RearrangeDeviceData(ListDevice,Wbase);

%%
% ==================================================
% Descriptor state space model
% ==================================================

% ### Power flow analysis
fprintf('Do the power flow analysis.\n')
switch Flag_PowerFlow
    case 1  % Gauss-Seidel 
        [PowerFlow,~,~,~,~,~,~,~] = PowerFlow_GS(ListBus,ListLine,Wbase);
    case 2  % Newton-Raphson
        [PowerFlow] = PowerFlow_NR(ListBus,ListLine,Wbase);
    otherwise
        error(['Error: Wrong setting for power flow algorithm.']);
end
ListPowerFlow = RearrangePowerFlow(PowerFlow);
% Move load flow to bus admittance matrix
[ListBus,ListLine,PowerFlow] = Load2SelfBranch(ListBus,ListLine,DeviceType,PowerFlow);
ListPowerFlow_ = RearrangePowerFlow(PowerFlow);

% ### Get the model of lines
fprintf('Get the descriptor-state-space model of network lines.\n')

[YbusObj,YbusDSS,~] = YbusCalcDSS(ListLine,Wbase);
[~,lsw] = size(YbusDSS.B);
ZbusObj = obj_SwitchInOut(YbusObj,lsw);
[ZbusStateStr,ZbusInputStr,ZbusOutputStr] = ZbusObj.ReadString(ZbusObj);

[st_YbusObj,st_YbusDSS,~] = YbusCalcDSS_Steady(ListLine,Wbase);
[~,lsw] = size(st_YbusDSS.B);
st_ZbusObj = obj_SwitchInOut(st_YbusObj,lsw);
[ZbusStateStr_st,ZbusInputStr_st,ZbusOutputStr_st] = st_ZbusObj.ReadString(st_ZbusObj);
% Notes:
% "st" means steady.

% ### Get the models of bus devices
fprintf('Get the descriptor-state-space model of bus devices.\n')
for i = 1:N_Device
    [GmObj_Cell{i},GmDSS_Cell{i},DevicePara{i},DeviceEqui{i},DeviceDiscreDamping{i},DeviceStateStr{i},DeviceInputStr{i},DeviceOutputStr{i}] = ...
        DeviceModel_Create('Type', DeviceType{i} ,'Flow',PowerFlow{i},'Para',Para{i},'Ts',Ts);
    
    % The following data is not used in the script, but will be used in
    % simulations, please do not delete.
    x_e{i} = DeviceEqui{i}{1};
    u_e{i} = DeviceEqui{i}{2};
    OtherInputs{i} = u_e{i}(3:end,:);
end

% ### Get the model of whole system
fprintf('Get the descriptor-state-space model of whole system.\n')
GmObj = DeviceModel_Link(GmObj_Cell);
[GsysObj,GsysDSS,Port_v,Port_i,Port_w,Port_T_m,Port_ang_r,Port_P_dc,Port_v_dc] = ...
    GmZbus_Connect(GmObj,ZbusObj);
[st_GsysObj,st_GsysDSS,st_Port_v,st_Port_i,st_Port_w,st_Port_T_m,st_Port_ang_r,st_Port_P_dc,st_Port_v_dc] = ...
    GmZbus_Connect(GmObj,st_ZbusObj);

% ### Chech if the system is proper
fprintf('Check if the whole system is proper:\n')
if isproper(GsysDSS)
    fprintf('Proper.\n');
    fprintf('Calculate the minimum realization of the system model for later use.\n')
    GminSS = minreal(GsysDSS);    
    % This "minreal" function only changes the element sequence of state
    % vectors, but does not change the element sequence of input and output
    % vectors.
    InverseOn = 0;
else
    error('Error: System is improper, which has more zeros than poles.')
end
if is_dss(GminSS)
    error(['Error: Minimum realization is in descriptor state space (dss) form.']);
end

% ### Output the System
fprintf('\n')
fprintf('==================================\n')
fprintf('Print the system\n')
fprintf('==================================\n')
fprintf('System object name: GsysObj\n')
fprintf('System name: GsysDSS\n')
fprintf('Minimum realization system name: GminSS\n')
if Enable_PrintOutput
    [SysStateString,SysInputString,SysOutputString] = GsysObj.ReadString(GsysObj);
    fprintf('Print power flow results: | bus | P | Q | V | angle | omega |')
    ListPowerFlow
    fprintf('Print system ports:\n')
    PrintSysString(N_Device,DeviceType,DeviceStateStr,DeviceInputStr,DeviceOutputStr,ZbusStateStr);
end
    
%%
% ==================================================
% Create Simulink Model
% ==================================================
fprintf('\n')
fprintf('=================================\n')
fprintf('Simulink Model\n')
fprintf('=================================\n')

if Enable_CreateSimulinkModel == 1
    
    fprintf('Create the simulink model aotumatically.\n')

    % Set the simulink model name
    Name_Model = 'mymodel_v1';

    % Close existing model with same name
    close_system(Name_Model,0);
    
    % Create the simulink model
    Main_Simulink(Name_Model,ListLine,DeviceType,ListAdvance,PowerFlow);
    fprintf('Get the simulink model successfully.\n')
    fprintf('Warning: for later use of the simulink model, please "save as" a different name.\n')

else
    fprintf('Warning: The auto creation of simulink model is disabled.\n')
end

%%
% ==================================================
% Plot
% ==================================================
    
fprintf('\n')
fprintf('==================================\n')
fprintf('Plot\n')
fprintf('==================================\n')

figure_n = 1000;

% Plot pole/zero map
fprintf('Calculate pole/zero.\n')
pole_sys = pole(GsysDSS)/2/pi;
st_pole_sys = pole(st_GsysDSS)/2/pi;
fprintf('Check if the system is stable:\n')
if isempty(find(real(pole_sys)>1e-9, 1))
    fprintf('Stable.\n');
else
    fprintf('Warning: Unstable.\n')
    Index = find(real(pole_sys)>1e-9);
    Unstable_Pole = pole_sys(Index)
end
if Enable_PlotPole
    fprintf('Plot pole/zero map.\n')
    figure_n = figure_n+1;
    figure(figure_n);
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    mtit('Global Pole Map: Actual');
    
	figure_n = figure_n+1;
    figure(figure_n);
    scatter(real(pole_sys),imag(pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    mtit('Zoomed Pole Map: Actual');
    axis([-250,50,-250,250]);
    
  	figure_n = figure_n+1;
    figure(figure_n);
    scatter(real(st_pole_sys),imag(st_pole_sys),'x','LineWidth',1.5); hold on; grid on;
    xlabel('Real Part (Hz)');
    ylabel('Imaginary Part (Hz)');
    mtit('Zoomed Pole Map: No Line EMT Dynamics');
    axis([-250,50,-250,250]);
else
    fprintf('Warning: The default plot of pole map is disabled.\n')
end

omega_p = logspace(-2,4,5e3)*2*pi;
omega_pn = [-flip(omega_p),omega_p];

% Plot admittance
if Enable_PlotAdmittance
    fprintf('Calculate complex-form admittance.\n')
    Tj = [1 1j;     % real to complex transform
          1 -1j];  
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            if InverseOn == 0
                Gr_ss{k} = GminSS(Port_i([2*k-1,2*k]),Port_v([2*k-1,2*k]));
            else          
                Gr_ss{k} = GminSS(Port_v([2*k-1,2*k]),Port_i([2*k-1,2*k]));
            end
            Gr_sym{k} = ss2sym(Gr_ss{k});
            Gr_c{k} = Tj*Gr_sym{k}*Tj^(-1);
        end
    end
    fprintf('Plot admittance.\n')
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if DeviceType{k} <= 50
            bodec(Gr_c{k}(1,1),1j*omega_pn,2*pi,'InverseOn',InverseOn,'PhaseOn',0); 
            CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)];
        end
    end
    legend(VecLegend);
    mtit('Bode Diagram: Admittance');
else
    fprintf('Warning: The default plot of admittance spectrum is disabled.\n')
end

% Plot w related
if Enable_PlotSwing
    fprintf('Find the w port relation.\n')
    for k = 1:N_Bus
        if floor(DeviceType{k}/10) == 0
            Gt_ss{k} = GminSS(Port_w(k),Port_T_m(k));
            Gt_sym{k} = -ss2sym(Gt_ss{k});
        elseif floor(DeviceType{k}/10) == 1
         	Gt_ss{k} = GminSS(Port_w(k),Port_ang_r(k));
            Gt_sym{k} = -ss2sym(Gt_ss{k});
        end
    end
    fprintf('Plot w port dynamics.\n')
 	figure_n = figure_n+1;
 	figure(figure_n);
    CountLegend = 0;
    VecLegend = {};
    for k = 1:N_Bus
        if (floor(DeviceType{k}/10) == 0) || (floor(DeviceType{k}/10) == 1)
            bodec(Gt_sym{k},1j*omega_pn,2*pi,'InverseOn',InverseOn,'PhaseOn',0);      
         	CountLegend = CountLegend + 1;
            VecLegend{CountLegend} = ['Bus',num2str(k)]; 
        end
    end
    legend(VecLegend);
    mtit('Bode Diagram: Swing');
    
else
    fprintf('Warning: The default plot of swing spectrum is disabled.\n');
end
    
%%
fprintf('\n')
fprintf('==================================\n')
fprintf('End: toolbox run successfully.\n')
fprintf('==================================\n')
   
