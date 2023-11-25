%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: shell.m
% By: Stephie Fried, David Lagakos
% Date: Winter 2022
% Purpose: Calibrates model for each country, for sensitivity exercise, and
% extension. Computes counterfactual experiments. Calculates tables and
% figures. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% Run the programs from the outermost folder in the replication package
% (click "add to path" after hitting run)
reppath = pwd; %Outermost folder in replication pacakage
bpath = strcat(reppath, '/matlab_programs'); %path to folder with matlab programs
tablePath = strcat(reppath, '/tables'); %path to folder where tables are stored
figpath = strcat(reppath, '/figures'); %path to folder where figures are stored
datapath = strcat(reppath, '/data/derived'); %path to folder where figures are stored

cd(bpath);

%Save paths
save ('paths', 'bpath', 'tablePath', 'figpath', 'reppath', 'datapath');

%Calibrate model
calibrate_shell

%Calculate steady states, decompositions, and weak links and tax
%experiments. 
steady_states

%Produce figures
figures

%Produce tables
tables
toc