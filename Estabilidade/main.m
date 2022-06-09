clear all
close all
clc
addpath('../aircraft-design-tool-main')

%% Read Excel
% Excel Name
filename = 'Centros_Geometricos.xlsx';

% get parts names
xlRange_comp = 'A2:A20';
[~,comp] = xlsread(filename,xlRange_comp);

% get parts json names
xlRange_jname = 'E2:E20';
[~,jname] = xlsread(filename,xlRange_jname);

% get coordinates
xlRange = 'B2:D20';
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
      0 -1  0];
rot_coord = R*[data_X_corr';data_Y_corr';data_Z_corr'];
X = rot_coord(1,:);
Y = rot_coord(2,:);
Z = rot_coord(3,:);

% Plot geometric Centers of components
plot_centers(X,Y,Z,comp)

clearvars -except X Y Z comp jname

%% Read Json file
% Get data from json and algorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('rescue.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

% Get mass for each component
comp_mass(1) = 0;
for i=2:length(jname)
    clear elem id
    if strcmp(jname{i}, 'Tail')
        jname_t1 = 'Vertical Tail';
        jname_t2 = 'Horizontal Tail';
        [elem1, id1] = find_by_name(data.vehicle.components, jname_t1);
        [elem2, id2] = find_by_name(data.vehicle.components, jname_t2);
        comp_mass(i) = elem1.mass + elem2.mass;
        clear jname_t1 jname_t2 elem1 elem2 id1 id2
    else
        [elem, id] = find_by_name(data.vehicle.components, jname{i});
        comp_mass(i) = elem.mass;
    end
end

%% Compute Center of Mass
% sum of masses of all components
Total_mass = sum(comp_mass(2:end));

% calculate Center of Mass
X_CG = sum(X(2:end).*comp_mass(2:end)) / Total_mass;
Y_CG = sum(Y(2:end).*comp_mass(2:end)) / Total_mass;
Z_CG = sum(Z(2:end).*comp_mass(2:end)) / Total_mass;

% Plot geometric Centers of components and Center of Mass
X_aux = X;
Y_aux = Y;
Z_aux = Z;
comp_aux = comp;
X_aux(end+1) = X_CG;
Y_aux(end+1) = Y_CG;
Z_aux(end+1) = Z_CG;
comp_aux{end+1} = 'CG';
plot_centers(X_aux,Y_aux,Z_aux,comp_aux)
