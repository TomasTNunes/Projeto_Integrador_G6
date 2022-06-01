clear all
close all
clc

h_cr = 1500;        % altitude de cruzeiro (m)
V_cr = 80;          % velocidade de cruzeiro (m/s)
d = 2;              % diâmetro máximo da fuselagem (m)
l = 8;              % comprimento da fuselagem (m)
S = 19.13;          % área da asa (m^2)
gamma = 1.4;        % coeficiente adiabático do ar
R = 287.052874;     % constante dos gases

T = 288.15-0.0065*h_cr;                         % temperatura (K)
rho = 1.225*(1-0.0065*h_cr/288.15)^4.25588;     % densidade (kg/m^3)
miu = 1.458*10^-6*T^(3/2)/(T+110.4);            % viscosidade (kg/ms)
M = sqrt(gamma*R*T);                            % número de Mach

q = 1/2*rho*V_cr^2;     % pressão dinâmica
Re = rho*l*V_cr/miu;    % Reynolds

A_side = 12335396.35*10^(-6);       % estimativa da área vista de lado
A_top = 14347845.53*10^(-6);        % estimativa da área vista de cima
S_wet = 1.7*(A_side+A_top);         % wetted area da fuselagem

form = 1 + 60*(l/d)^3 + (d/l)/400;                      % form factor
Q = 1.5;                    % interference factor igual ao do exemplo

Cf = 0.455/(((log10(Re))^2.58)*(1+0.144*M^2)^0.65); % coeficiente de fricção turbulento

Ff = q*S_wet*Cf*form*Q;
Fw = 0;                     % não há wave drag

Cd0 = (Ff+Fw)/(q*S)   

%-------------------------------------------------------------------

MTOM = 2187.59/0.45359;       % maximum take-off mass/weight (lbs)
n = 3;                        % número de rodas do trem de aterragem principal

W_wheel = 0.9*MTOM/n;         % peso suportado por roda (lbs)

d_wheel = 0.0254*1.51*W_wheel^0.349       % diâmetro da roda (m)
width_wheel = 0.0254*0.715*W_wheel^0.312  % largura da roda (m)

% Assumir trem de aterragem do nariz 40% mais pequeno?



