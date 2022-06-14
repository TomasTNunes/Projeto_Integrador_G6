clear all
close all
clc
addpath('../aircraft-design-tool-main')

%% Read Excel
% Excel Name
filename = 'Centros_Geometricos.xlsx';

% get parts names
xlRange_comp = 'A2:A18';
[~,comp] = xlsread(filename,xlRange_comp);

% get parts json names
xlRange_jname = 'E2:E18';
[~,jname] = xlsread(filename,xlRange_jname);

% get program frame
xlRange_prog_frame = 'F2:F18';
[~,prog_frame] = xlsread(filename,xlRange_prog_frame);

% get coordinates
xlRange = 'B2:D18';
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
data_Z = data(:,3);

%% Data Processing
% Translate Frame to chosen origin
data_X_corr = data_X(strcmp(prog_frame,'SW')) - data_X(1);
data_X_corr = [data_X_corr(1);data_X(2);data_X(3);data_X_corr(2:end)];
data_Y_corr = data_Y(strcmp(prog_frame,'SW')) - data_Y(1);
data_Y_corr = [data_Y_corr(1);data_Y(2);data_Y(3);data_Y_corr(2:end)];
data_Z_corr = data_Z(strcmp(prog_frame,'SW')) - data_Z(1);
data_Z_corr = [data_Z_corr(1);data_Z(2);data_Z(3);data_Z_corr(2:end)];

% Rotate from SW frame to body frame 
R = [-1  0  0;
      0  0 -1;
      0 -1  0];
rot_coord = R*[data_X_corr(strcmp(prog_frame,'SW'))';data_Y_corr(strcmp(prog_frame,'SW'))';data_Z_corr(strcmp(prog_frame,'SW'))'];
X = rot_coord(1,:);
Y = rot_coord(2,:);
Z = rot_coord(3,:);

% Rotate from XFLR5 frame to body frame 
X = [X(1),-data_X_corr(2),-data_X_corr(3),X(2:end)];
Y = [Y(1),data_Y_corr(2),data_Y_corr(3),Y(2:end)];
Z = [Z(1),-data_Z_corr(2),-data_Z_corr(3),Z(2:end)];

% Plot geometric Centers of components
plot_centers(X,Y,Z,comp)

%clearvars -except X Y Z comp jname

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
num_fuelt = strcmp(jname, 'Fuel Tank');
num_fuelt = length(num_fuelt(num_fuelt==1));
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
    elseif strcmp(jname{i}, 'Fuel Tank')
        [elem, id] = find_by_name(data.vehicle.components, jname{i});
        comp_mass(i) = elem.mass / num_fuelt;
    elseif strcmp(jname{i}, 'Fuselage')
        [elem, id] = find_by_name(data.vehicle.components, jname{i});
        comp_mass(i) = elem.mass + data.vehicle.components{1, 1}.mass + ...
            data.vehicle.components{2, 1}.mass + data.vehicle.components{4, 1}.mass;
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

% XFLR-5 frame
pos_xflr5 = [-X' Y' -Z'];
