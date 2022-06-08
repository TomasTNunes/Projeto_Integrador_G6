clear all
close all
clc
addpath('../aircraft-design-tool-main')

%% Read Excel
% Excel Name
filename = 'Centros_Geometricos.xlsx';

% get parts names
xlRange_comp = 'A2:A16';
[~,comp] = xlsread(filename,xlRange_comp);
% get parts json names
xlRange_jname = 'E2:E16';
[~,jname] = xlsread(filename,xlRange_jname);
% get coordinates
xlRange = 'B2:D16';
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
data_Z = data(:,3);

%% Data Processing
% Translate Frame to chosen origin
data_X_corr = data_X - data_X(1);
data_Y_corr = data_Y - data_Y(1);
data_Z_corr = data_Z - data_Z(1);

% Rotate from SW frame to body frame
R = [-1  0  0;
      0  0 -1;
      0 -1  0;];
rot_coord = R*[data_X_corr';data_Y_corr';data_Z_corr'];
X = rot_coord(1,:);
Y = rot_coord(2,:);
Z = rot_coord(3,:);

% Plot geometric Centers of components
plot_geo_centers(X,Y,Z,comp)

clearvars -except X Y Z comp jname

%% Read Json file
% Get data from json and algorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('transport.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

