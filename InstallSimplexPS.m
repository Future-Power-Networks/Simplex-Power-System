% This file installs the SimplexPowerSystem toolbox on the users' PC. Users
% do not need to read or understand the codes in this file. Please just run
% it directly and run it only ONCE.

% Author(s): Yitong Li, Yunjie Gu

% Notes for developers:
% The simulink library may need to be updated programmatically here, when
% adding Simplex Library.

%%
clear all
clc
close all

%%
% Change the current folder
fprintf('Changing the current folder to the toolbox folder...\n')
mfile_name = mfilename('fullpath');
[RootPath,~,~]  = fileparts(mfile_name);
cd(RootPath);

% Check the matlab version
fprintf('Checking the Matlab version...\n')
MatlabVersion = version('-release');
MatlabVersion = MatlabVersion(1:(end-1));
MatlabVersionYear = MatlabVersion;
if str2double(MatlabVersionYear)<2015
    error(['Error: Please use Matlab version 2015a or later!']);
end

% Check if a previous version of SimplexPS has been installed
fprintf('Checking if the SimplexPS has been installed before...\n')
if exist('SimplexPS')~=0
    error(['Error: SimplexPowerSystem has been installed on this PC/laptop before. Please unstall the old version of SimplexPowerSystem first!']);
end

% Add folder to path
fprintf('Installing SimplexPowerSystem...\n')
addpath(genpath(pwd));  % Add path
savepath;

% Convert the toolbox lib to the required version
fprintf('Converting the toolbox library to the required Matlab version, please wait a second...\n')
warning('off','all')    % Turn off all warnings
load_system('SimplexPS_2015a.slx');
save_system('SimplexPS_2015a.slx','Library/SimplexPS.slx');
close_system('SimplexPS.slx');
warning('on','all')     % Turn on all warnings
clc

%%
% Installation is completed
DlgTitle = 'Congratulations!';
DlgQuestion = 'SimplexPowerSystem is installed successfully! Do you want to run "UserMain.m" now to use the toolbox?';
choice = questdlg(DlgQuestion,DlgTitle,'Yes','No','Yes');

%%
% Notes: 
% adding path should be improved.

if strcmp(choice,'Yes')
    open('UserMain.m');
    run('UserMain.m');
else
 	msgbox('The installation of SimplexPowerSystem is completed. Please run "UserMain.m" later in the root path for using the toolbox!');
end