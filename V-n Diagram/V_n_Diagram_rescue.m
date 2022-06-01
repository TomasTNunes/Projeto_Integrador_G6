clear all
close all
clc
addpath('../aircraft-design-tool-main')

% Get data from json and anlgorithm
global constants;
constants.g = 9.81; % m/s^2
data = load_project('rescue.json');
data.mission = build_mission(data.mission);
data.vehicle = build_vehicle(data.mission, data.vehicle);
data.vehicle = aero_analysis(data.mission, data.vehicle);
[data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);

% Parameters
MTOW = data.vehicle.mass;
W = MTOW * constants.g;
S = data.vehicle.components{6, 1}.aspect_ratio * (data.vehicle.components{6, 1}.mean_chord)^2;
CLmax = data.vehicle.components{6, 1}.airfoil.cl_max  ;
CLalpha = data.vehicle.components{6, 1}.airfoil.lift_slope_coefficient;
c = data.vehicle.components{6, 1}.mean_chord;
Vc = data.mission.segments{6, 1}.velocity;
rho = data.mission.segments{6, 1}.density; % h=700m
n_max = 3;
n_min = -1;

% V-n relation
n_func = @(v) v.^2.*0.5*rho*S*CLmax/W;
v_func = @(n) (n.*2*W/(S*rho*CLmax)).^0.5;

% Velocity for limits
Va = v_func(n_max);
Vg = v_func(abs(n_min));
Vd = 1.5*Vc;

% Positive Stall Limit
V_PSL = linspace(0,Va,100);
n_PSL = n_func(V_PSL);

% Negative Stall Limit
V_NSL = linspace(0,Vg,100);
n_NSL = -n_func(V_NSL);

% Upper Structural Limit
V_USL = linspace(Va,Vd,100);
n_USL = zeros(length(V_USL))+n_max;

% Lower Structural Limit
V_LSL = linspace(Vg,Vc,100);
n_LSL = zeros(length(V_LSL))+n_min;

% Aeroelastic Limit
m_AL = (0-n_min)/(Vd-Vc);
b_AL = -m_AL*Vd;
n_AL_func = @(v) m_AL.*v+b_AL;
V_AL = linspace(Vc,Vd,100);
n_AL = n_AL_func(V_AL);

% Dynamic pressure Limit
V_DPL = [Vd Vd];
n_DPL = [n_max 0];

% Plot Basic V-n Diagram
figure(1)
plot(V_PSL,n_PSL,'-r','LineWidth',1.5)
hold on
plot(V_NSL,n_NSL,'-r','LineWidth',1.5)
plot(V_USL,n_USL,'-r','LineWidth',1.5)
plot(V_LSL,n_LSL,'-r','LineWidth',1.5)
plot(V_AL,n_AL,'-r','LineWidth',1.5)
plot(V_DPL,n_DPL,'-r','LineWidth',1.5)

% delta n due to gust
dn_func = @(u,v) 0.5*rho*u*v*CLalpha*S/W;

% Mass ratio and constant K
miu = 2*W/(S*rho*constants.g*c*CLalpha);
K = 0.88*miu/(5.3+miu); % subsonic (M<1)

% Gust velocities
U1 = K*66*0.3048; % for max maneuvering speed (Va)
U2 = K*50*0.3048; % for cruise speed (Vc)
U3 = K*25*0.3048; % for dive speed (Vd)

% Calculate respectives delta n and n
dn1 =  dn_func(U1,Va);
dn2 =  dn_func(U2,Vc);
dn3 =  dn_func(U3,Vd);

n1_p = 1+dn1;
n1_n = 1-dn1;
n2_p = 1+dn2;
n2_n = 1-dn2;
n3_p = 1+dn3;
n3_n = 1-dn3;

% Plot V-n Diagram with gust
V_wG = [0 Va Vc Vd Vd Vc Va 0];
n_wG = [1 n1_p n2_p n3_p n3_n n2_n n1_n 1];
plot(V_wG,n_wG,'-b','LineWidth',1.5)

% Determine Combined V-n Diagram
V_n1 = v_func(1); % V(n=1)

V_1 = linspace(V_n1,Va,100);
n_1 = n_func(V_1);

m_2 = (n2_p-n_max)/(Vc-Va);
b_2 = n_max-m_2*Va;
n_2_func = @(v) m_2.*v +b_2;
V_2 = linspace(Va,Vc,100);
n_2 = n_2_func(V_2);

m_3 = (Vc-Vd)/(n2_p-n3_p);
b_3 = Vd-m_3*n3_p;
V_3_func = @(n) m_3.*n +b_3;
V_3 = [Vc V_3_func(n_max)];
n_3 = [n2_p n_max];

V_4 = [V_3_func(n_max) Vd];
n_4 = [n_max n_max];

V_5 = [Vd Vd];
n_5 = [n_max n3_n];

V_6 = [Vd Vc];
n_6 = [n3_n n2_n];

V_7 = [Vc Va];
n_7 = [n2_n n1_n];

m_8 = (Va-0)/(n1_n-1);
b_8 = 0-m_8*1;
V_8_func = @(n) m_8.*n +b_8;
V_8 = [Va V_8_func(n_min)];
n_8 = [n1_n n_min];

V_9 = [V_8_func(n_min) V_n1];
n_9 = [n_min n_min];

V_10 = [V_n1 V_n1];
n_10 = [n_min 1];

% Duplicate figure
a1 = gca;
f2 = figure(2);
a2 = copyobj(a1,f2);
f=get(gca,'Children');
legend([f(1),f(2)],'Basic V-n Diagram','V-n Diagram with gust','Location','Northwest')
title('V-n Diagram')
xlabel('V [m/s]')
ylabel('n')
xlim([0 Vd+10])
clear f

% Plot Combined V-n Diagram
figure(1)
plot(V_1,n_1,'-g','LineWidth',1.5)
plot(V_2,n_2,'-g','LineWidth',1.5)
plot(V_3,n_3,'-g','LineWidth',1.5)
plot(V_4,n_4,'-g','LineWidth',1.5)
plot(V_5,n_5,'-g','LineWidth',1.5)
plot(V_6,n_6,'-g','LineWidth',1.5)
plot(V_7,n_7,'-g','LineWidth',1.5)
plot(V_8,n_8,'-g','LineWidth',1.5)
plot(V_9,n_9,'-g','LineWidth',1.5)
plot(V_10,n_10,'-g','LineWidth',1.5)

% Plot Area in Combined V-n Diagram
patch('XData',[V_1,V_2,V_3,V_4,V_5,V_6,V_7,V_8,V_9,V_10],'YData',[n_1,n_2,n_3,n_4,n_5,n_6,n_7,n_8,n_9,n_10],'FaceAlpha',0.2,'FaceColor','green','LineStyle','none')

% Plot Info
f=get(gca,'Children');
legend([f(13),f(12),f(2)],'Basic V-n Diagram','V-n Diagram with gust','Combined V-n Diagram','Location','Northwest')
title('V-n Diagram')
xlabel('V [m/s]')
ylabel('n')
xlim([0 Vd+10])