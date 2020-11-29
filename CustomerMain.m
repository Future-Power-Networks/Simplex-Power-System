% To customer:
% Please use and only use this file to run toolbox.

%% Tips

% Please ensure that the toolbox is installed (by using "Install.m") the
% first time.

% For changing default data, please use "CustomerData.xlsx". (More examples
% can be found in the "Examples" folder)

% The toolbox defaultly prints the results in Matlab command window, saves
% the results into Matlab workspace, and prints figures.

%% Run toolbox
SimplexPS.Toolbox.Main(); 	% This function runs toolbox.

%% Custormer available results
GsysDSS;            % Whole-system model (descriptor state space form)
                    % Note: The elements of state, input, and output
                    % vectors are printed in the command window.
                    
GminSS;             % Whole-system model (state space form)
                    % Note: This model is the minimum realization of
                    % GsysDSS, which keeps the same input and output as
                    % GsysDSS, but reduces the order of state.
                    
ListPowerFlow;      % Power flow result
                    % Note: The result is in the form of
                    % | bus | P | Q | V | angle | omega |

ListPowerFlow_;     % Power flow result only for active device by combing 
                    % the load into the nodal admittance matrix.
                    
pole_sys;           % Whole-system poles, or equivalently eigenvalues.

%mymodel_v1;         % This is the simulink model generated automatically 
                    % based on "CustomerData.xlsx".

%% Customer function and plot
% Custormer can write their functions here to further deal with the data
% mentioned above.