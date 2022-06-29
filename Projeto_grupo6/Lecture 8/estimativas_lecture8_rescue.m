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

% Parameters
h_cr = data.mission.segments{6, 1}.altitude;        % altitude de cruzeiro (m)
V_cr = data.mission.segments{6, 1}.velocity;        % velocidade de cruzeiro (m/s)
d = data.vehicle.components{5, 1}.diameter;         % diâmetro máximo da fuselagem (m)
l = data.vehicle.components{5, 1}.length;           % comprimento da fuselagem (m)
S = data.vehicle.components{6, 1}.aspect_ratio*(data.vehicle.components{6, 1}.mean_chord)^2;;     % área da asa (m^2)

T = data.mission.segments{6, 1}.temperature;    % temperatura (K)
rho = data.mission.segments{6, 1}.density;     % densidade (kg/m^3)
a = data.mission.segments{6, 1}.speed_sound;   % Velocidade do Som
miu = 1.458*10^-6*T^(3/2)/(T+110.4);           % viscosidade (kg/ms)

q = 1/2*rho*V_cr^2;     % pressão dinâmica
Re = rho*l*V_cr/miu;    % Reynolds

A_side = 12335396.35*10^(-6);       % estimativa da área vista de lado
A_top = 14347845.53*10^(-6);        % estimativa da área vista de cima
S_wet = 1.7*(A_side+A_top);         % wetted area da fuselagem

% form = 1 + 60*(l/d)^3 + (d/l)/400;            % fórmula dos slides errada
form = 1 + 60*(d/l)^3 + (l/d)/400;              % fórmula do exemplo
Q = data.vehicle.components{5, 1}.interf_factor;         % interference factor igual ao do exemplo      

Cf = 0.455/(((log10(Re))^2.58)*(1+0.144*(V_cr/a)^2)^0.65); % coeficiente de fricção turbulento

Ff = q*S_wet*Cf*form*Q;
Fw = 0;                     % não há wave drag

Cd0_fuselage = (Ff+Fw)/(q*S)   

%-------------------------------------------------------------------

MTOM = data.vehicle.mass/0.45359;    % maximum take-off mass/weight (lbs)
n = 3;                        % número de rodas do trem de aterragem principal

W_wheel = 0.9*MTOM/n;         % peso suportado por roda (lbs)

d_wheel = 0.0254*1.51*W_wheel^0.349       % diâmetro da roda (m)
width_wheel = 0.0254*0.715*W_wheel^0.312  % largura da roda (m)

% Assumir trem de aterragem do nariz 40% mais pequeno?