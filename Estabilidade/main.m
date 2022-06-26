clear all
close all
clc
addpath('../aircraft-design-tool-main')

%% Read Excel
% Excel Name
% filename = 'Centros_Geometricos.xlsx';
filename = 'Centros_Geometricos_estV.xlsx';

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
%plot_centers(X,Y,Z,comp)

clearvars -except X Y Z comp jname

%% Read Json file
% Get data from json and algorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('rescue_estV.json'); % rescue/rescue_estV
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
    elseif strcmp(jname{i}, 'Fuel Tank')
        [elem, id] = find_by_name(data.vehicle.components, jname{i});
        comp_mass(i) = elem.mass * (1 + elem.reserve);
    elseif strcmp(jname{i}, 'Battery')
        [elem, id] = find_by_name(data.vehicle.components, jname{i});
        comp_mass(i) = elem.mass * (1 + elem.reserve);
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
X_CG = sum(X(2:end).*comp_mass(2:end)) / Total_mass
Y_CG = sum(Y(2:end).*comp_mass(2:end)) / Total_mass
Z_CG = sum(Z(2:end).*comp_mass(2:end)) / Total_mass

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

%% Vertical Stability
mission = data.mission;
vehicle = data.vehicle;
energy = data.energy;

% Get parameters
MTOM = vehicle.mass;
Vtip = vehicle.components{12, 1}.tip_velocity; % Rotor = Propeller
r_rot = vehicle.components{12, 1}.radius; % Rotor = Propeller
ki = vehicle.components{12, 1}.induced_power_factor; % Rotor induced power factor
cd0 = vehicle.components{12, 1}.base_drag_coefficient; % Rotor base drag coefficient
s = vehicle.components{12, 1}.rotor_solidity; % Rotor solidity (-)
rho_vtol = mission.segments{3, 1}.density; % h=150m
NR = vehicle.components{12, 1}.number; % Number of Rotors
[~, DL, ~, ~] = calculates_loadings(mission, vehicle, energy); % (N/m^2) Disk loading
Vy = mission.segments{2, 1}.velocity; % Vertical Climb velocity

% Required Thrust for Hover and Vertcical Climb
Tr = MTOM*constants.g; % Thrust in hover (N)

% Compute Rotors Thrust so that X_CG = X of the center of Forces from rotor
xa = pos_xflr5(4,1); % x position of wing rotors
xc = pos_xflr5(6,1); % x position of tail rotors
Tc = @(Ta) (Tr-2*Ta)/2; % Ta -> Thrust in one wing rotor & Tc -> Thrust in one tail rotor
M = @(Ta) (xc+X_CG)*Tc(Ta)*2 - (-X_CG-xa)*Ta*2; % Torque in CG
options = optimset('Display','off');
Ta = fsolve(M,Tr/4,options) % Solve M=0 in order to Ta (Thrust in one wing rotor) (N)
Tc = Tc(Ta) %Thrust in one tail rotor (N)

% Verify if these Thrusts are within Power Limits of eletric motor
Phra = Ta*(ki*sqrt(DL/(2*rho_vtol)) + (rho_vtol*Vtip*Vtip*Vtip/DL)*(s*cd0/8)) % Hover power in one wing rotor (W)
Pclra = Ta*(Vy - ki*Vy/2 + ki*0.5*sqrt(Vy*Vy + 2*DL/rho_vtol) + ((rho_vtol*Vtip*Vtip*Vtip)/DL)*(s*cd0/8)) % Power in vertical climb in one wing rotor (W)
Phrc = Tc*(ki*sqrt(DL/(2*rho_vtol)) + (rho_vtol*Vtip*Vtip*Vtip/DL)*(s*cd0/8)) % Hover power in one tail rotor (W)
Pclrc = Tc*(Vy - ki*Vy/2 + ki*0.5*sqrt(Vy*Vy + 2*DL/rho_vtol) + ((rho_vtol*Vtip*Vtip*Vtip)/DL)*(s*cd0/8)) % Power in vertical climb in one tail rotor (W)
