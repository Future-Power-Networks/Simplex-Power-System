% To users:
% Please use and ONLY use this file to run toolbox.

%% Tips

% Please ensure that the toolbox is installed， by running
% "InstallSimplexPS.m" the first time.

% The toolbox defaultly prints the results in Command Window, saves the
% results into Workspace, and plots figures.

% For changing default user data, please change "UserData.xlsx". More
% examples can be found in the "Examples" folder.

%% Clear matlab
clear all;  % Clear matlab workspace
clc;        % Clear matlab command window
close all;  % Close all figures, etc

%% Set user data
% Default
Name_Netlist = 'UserData.xlsx';
% "UserData.xlsx" defaultly contains the data of a 4-bus
% generator-inverter-composite power system. Please feel free to change it.

% Other standard example systems ("Examples" folder):
% Name_Netlist = 'IEEE_14Bus';
% Name_Netlist = 'IEEE_30Bus';      % To be added
% Name_Netlist = 'IEEE_57Bus';     	% To be added
% Name_Netlist = 'IEEE_118Bus';     % To be added

% Examples used in papers (subfolers in "Examples" folder):
% Name_Netlist = 'PortCoupling_SingleSGInfiniteBus.xlsx';  	% Single-generator-infinite-bus system
% Name_Netlist = 'PortCoupling_SingleVSIInfiniteBus.xlsx';	% Single-inverter-infinite-bus system
% Name_Netlist = 'PortCoupling_IEEE14Bus.xlsx';            	% 14 bus system
% Name_Netlist = 'Duality_GridFormingInverter';
% Name_Netlist = 'Duality_GridFollowingInverter';
% Name_Netlist = 'Duality_14Bus';

%% Run toolbox
SimplexPS.Toolbox.Main();

%% Results available to users
% GsysDSS;          % Whole-system port model (descriptor state space
                    % form). Notes: The elements of state, input, and
                    % output vectors are printed in the command window.
                    %
                    % A quick introduction of DSS modeling method:
                    % https://uk.mathworks.com/help/simulink/slref/descriptorstatespace.html
                    
% GminSS;           % Whole-system port model (state space form).
                    % Notes: This model is the minimum realization of
                    % GsysDSS, which keeps the same input and output as
                    % GsysDSS, but reduces the order of state.

% YsysDSS;          % Whole-system admittance model (descriptor state space
                    % form). Notes: This model is derived from GsysDSS by
                    % keeping the voltage and current ports only and
                    % removing other inputs and outputs.
                    
% ListPowerFlow;    % Power flow result
                    % Notes: The result is in the form of
                    % | bus | P | Q | V | angle | omega |

% ListPowerFlow_;   % Power flow result only for active device by combing 
                    % the load into the nodal admittance matrix.
                    
% pole_sys;         % Whole-system poles, or equivalently eigenvalues.

% mymodel_v1;       % This is the simulink model generated automatically 
                    % based on the user data.

%% User function and plot
% Users can write their own functions here to further deal with the data
% mentioned above.