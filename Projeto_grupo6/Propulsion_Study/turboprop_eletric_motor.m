clear all
close all
clc
addpath('../aircraft-design-tool-main/')

% Get data from json and algorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('rescue.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

mission = data.mission;
vehicle = data.vehicle;
energy = data.energy;

% Parameters [rescue]
MTOM = vehicle.mass;
[wl_design, DL, ~, ~] = calculates_loadings(mission, vehicle, energy);
S = MTOM * constants.g / wl_design;
Vc = mission.segments{6, 1}.velocity;
rho_cr = mission.segments{6, 1}.density; % h=700m
a_cr = mission.segments{6, 1}.speed_sound; % h=700m
cdo_cr = data.vehicle.segments{6, 1}.base_drag_coefficient; % cruise
e = 0.8;
AR = vehicle.components{6, 1}.aspect_ratio;
NP = vehicle.components{13, 1}.number; % Number of Propellers
ef_prop = vehicle.components{13, 1}.efficiency;
Vtip = vehicle.components{12, 1}.tip_velocity; % Rotor = Propeller
r_rot = vehicle.components{12, 1}.radius; % Rotor = Propeller
ki = vehicle.components{12, 1}.induced_power_factor; % Rotor induced power factor
cd0 = vehicle.components{12, 1}.base_drag_coefficient; % Rotor base drag coefficient
s = vehicle.components{12, 1}.rotor_solidity; % Rotor solidity (-)
rho_vtol = mission.segments{3, 1}.density; % h=150m
NR = vehicle.components{12, 1}.number; % Number of Rotors


% Turboprop
CD_cr = cdo_cr + 1/(pi*AR*e)*(MTOM*constants.g/(0.5*rho_cr*Vc^2*S));
D_cr = 0.5*rho_cr*Vc^2*S*CD_cr;
Tcr = D_cr;
Pcr = Tcr*Vc/ef_prop % more power need due to design point
% 2350 --> Max RPM of propeller
% RPM for cruise
n_cr = sqrt((Vtip^2 - Vc^2))/(pi*2*r_rot) *60 % RPM


% Hover (Out of Ground Effect)
Th = MTOM*constants.g; % Thrust in hover (N)
Thr = Th/NR; % Thrust in hover per rotor (N) - Assumption: all rotors contribute equally to the total thrust
Phr = Thr*(ki*sqrt(DL/(2*rho_vtol)) + (rho_vtol*Vtip*Vtip*Vtip/DL)*(s*cd0/8)) % Hover power per rotor (W)
Ph = NR*Phr; % Hover power (W)
% 2350 --> Max RPM of rotor
% RPM for hover
n_h = Vtip/(pi*2*r_rot) *60 % RPM


% Vertical Climb
Vy = mission.segments{2, 1}.velocity; % Vertical Climb velocity
Tclr = Thr;  % Thrust in vertical climb per rotor (N)- Assumption: the same thrust as in hover
Pclr = Tclr*(Vy - ki*Vy/2 + ki*0.5*sqrt(Vy*Vy + 2*DL/rho_vtol) + ((rho_vtol*Vtip*Vtip*Vtip)/DL)*(s*cd0/8)) % Power in vertical climb per rotor (W)
Pcl = NR*Pclr;  % Power in vertical climb (W)
