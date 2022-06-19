clear all
close all
clc
addpath('../aircraft-design-tool-main/')

% Get data from json
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

% Noise estimation
m = 1; % (-) Harmonic number
NB = vehicle.components{12, 1}.number_blades; % (-) Number of blades 
c = 0.35; % (m) Blade's mean chord
t = 0.12*c; % (m) Blade's maximum thickness (Assumption)
a = mission.segments{6, 1}.speed_sound; % (m/s) Speed of Sound
Vtip = vehicle.components{12, 1}.tip_velocity; % (m/s) Tip Velocity
r = vehicle.components{12, 1}.radius; % (m) Rotor Radius
Omega = Vtip/r; % (rad/s)

NR = vehicle.components{12, 1}.number; % Number of Rotors
Th = vehicle.mass * constants.g; % Thrust in hover (N)
Thr = Th / NR; % Thrust in hover per rotor (N) - Assumption: all rotors contribute equally to the total thrust
ki = vehicle.components{12, 1}.induced_power_factor; % Rotor induced power factor
[~, DL, fpl_design, ~] = calculates_loadings(mission, vehicle, energy); % (N/m^2) Disk loading
rho = mission.segments{6, 1}.density; % Air Density (kg/m3) - SL, ISA + 20ºC conditions [Cruise, mas a altura do cruise e baixa 700m]
s = vehicle.components{12, 1}.rotor_solidity; % Rotor solidity (-)
cd0 = vehicle.components{12, 1}.base_drag_coefficient; % Rotor base drag coefficient
Phr = Thr*(ki*sqrt(DL/(2*rho)) + (rho*Vtip*Vtip*Vtip/DL)*(s*cd0/8)); % Hover power per rotor (W)
Torque = Phr*r/Vtip; % (N.m)
Cl = 6*Thr/(rho*Vtip*Vtip*NB*c*r);  % (-) mean lift coefficient
alt = 500*0.3048; % (m) - vertical distance of the observer from the source, 500 ft
pref = 2.0e-5; % (Pa) Reference pressure
St = 0.28; % (-) Strouhal number
nsteps = 50;

for i=1:(nsteps+1)
    aa = 89 - (i-1)*80/nsteps;
    ydist(i) = alt/sin(aa*pi()/180);
    DS = ydist(i);
    tt(i) = (pi()/2) + aa*(pi()/180);
    theta = tt(i);

    % Rotational Noise
    Re = 0.8*r;
    JmB = besselj(m*NB,(m*NB*Omega*Re*sin(theta)/a));
    pmL = (m*NB*Omega/(2*sqrt(2)*pi()*a*DS))*(Thr*cos(theta) - Torque*a/(Omega*Re*Re))*JmB;
    pmT = -(rho*((m*NB*Omega)^2)*NB/(3*sqrt(2)*pi()*DS))*c*t*Re*JmB;
    pd_rot = NR*((pmL*pmL + pmT*pmT)/(pref*pref));

    % Vortex Noise
    AA = pi()*r*r; % Disk area
    % 1st Step - Find Peak Frequency and the Overall SPL
    alpha = Cl/(2*pi());
    hh = t*cos(alpha) + c*sin(alpha);
    fpeak = (0.7*Vtip)*St/hh;
    k2 = 1.206*10^-2;   % (s^3/ft^3) Calibrated for helicopter
    k2 = k2/(.3048^3);  % (s^3/m^3)
    SPLo = 20*log10(k2*(Vtip/(rho*DS))*sqrt(NR*Thr*(Thr/AA)/s));
    % 2nd Step - Calculate SPL for each of spectrum frequencies
    freq = [fpeak/2; fpeak; 2*fpeak; 4*fpeak; 8*fpeak; 16*fpeak];
    SPL_w = [7.92; 4.17; 8.33; 8.75; 12.92; 13.33];
    SPL_ctr = 0;
    for ii = 1:5
      fr1 = freq(ii)/fpeak;
      SPL1 = SPLo - SPL_w(ii);
      fr2 = freq(ii+1)/fpeak;
      SPL2 = SPLo - SPL_w(ii+1);
      C1 = (SPL2 - SPL1)/(log10(fr2) - log10(fr1));
      C2 = SPL2 - C1*log10(fr2);
      int_val = (10^(C2/10))*((fr2^(C1/10 + 1))/(C1/10 + 1)) - (10^(C2/10))*((fr1^(C1/10 + 1))/(C1/10 + 1));
      SPL_ctr = SPL_ctr + int_val;
    end
    pd_vortex = SPL_ctr;

    % Total Noise
    SPL(i) = 10*log10(pd_rot + pd_vortex);
end
SPL_max = max(SPL) % (dB)
