clear all
close all
clc

global constants;
constants.g = 9.81; % m/s^2

%% Raio Rotores, Nº Rotores
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(0.7,1,15);
Y_values = [2 4 6 8];
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{12, 1}.radius = X_values(i);
        data.vehicle.components{12, 1}.number = Y_values(ii);
        data.vehicle.components{16, 1}.number = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Raio dos Rotores [m]","Nº de Rotores","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','A1');

%% h vertical climb, Velocidade de subida vertical
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(50,250,15);
Y_values = linspace(1,6,15);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.mission.segments{2, 1}.altitude(2) = X_values(i);
        data.mission.segments{3, 1}.altitude = X_values(i);
        data.mission.segments{4, 1}.altitude = X_values(i);
        data.mission.segments{5, 1}.altitude(1) = X_values(i);
        data.mission.segments{2, 1}.velocity = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Altitude Subida Vertical [m]","Velocidade Subida Vertical [m/s]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','I1');

%% Densidade da bateria, Massa bateria
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(100000,2000000,100);
j = 1;

for i=1:length(X_values)
    data.vehicle.components{10, 1}.specific_energy = X_values(i);
    
    %adt
    data.mission = build_mission(data.mission);
    data.vehicle = build_vehicle(data.mission, data.vehicle);
    data.vehicle = aero_analysis(data.mission, data.vehicle);
    [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
    X(j,1) = X_values(i);
    Y(j,1) = data.vehicle.components{10, 1}.mass;
    j = j + 1;
end

%table
varNames = ["Densidade da bateria [J/kg]","Massa da bateria [kg]"];
tab = table(X,Y,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','Q1');

%% Aspect Ratio, Mean Chord
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(5.5,10,15);
Y_values = linspace(1.3,1.9,15);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{6, 1}.aspect_ratio = X_values(i);
        data.vehicle.components{6, 1}.mean_chord = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Aspect Ratio","Corda media [m]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','T1');

%% Crew pax payload, massa fuselagem
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(20,200,20);
Y_values = linspace(200,1000,20);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{4, 1}.mass = X_values(i);
        data.vehicle.components{5, 1}.mass = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Payload [kg]","Massa fuselagem [kg]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','AB1');

%% massa fuselagem, range
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(200,1000,20);
Y_values = linspace(50,500,20);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{5, 1}.mass = X_values(i);
        data.mission.segments{6, 1}.range = Y_values(ii)*1000;
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Massa fuselagem [kg]","Cruise range [km]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','AJ1');

%% bSFC cruise
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(2,20,100);
j = 1;

for i=1:length(X_values)
    data.energy.networks{2, 1}.layout{2, 1}.brake_specific_fuel_consumption = X_values(i)*10^(-8);
            
    %adt
    data.mission = build_mission(data.mission);
    data.vehicle = build_vehicle(data.mission, data.vehicle);
    data.vehicle = aero_analysis(data.mission, data.vehicle);
    [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
    [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
    X(j,1) = X_values(i);
    MTOW(j,1) = data.vehicle.mass;
    j = j + 1;
end

%table
varNames = ["bSFC cruise *10^-^8 [kg/Ws]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','AR1');

%% range cruise, v=cte
clearvars -except constants
data = load_project('transport.json');
X_values = 80;
Y_values = linspace(50,500,100);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.mission.segments{4, 1}.velocity(2) = X_values;
        data.mission.segments{5, 1}.velocity = X_values;
        data.mission.segments{6, 1}.velocity = X_values;
        data.mission.segments{6, 1}.range = Y_values(ii)*1000;
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values;
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Velocidade cruise [m/s]","Cruise range [km]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','AY1');

%% Nº propellers, Potencia turboprop
clearvars -except constants
data = load_project('transport.json');
X_values = [1 2 3];
Y_values = linspace(400000,1500000,15);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{13, 1}.number = X_values(i);
        data.vehicle.components{9, 1}.max_power = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Nº propellers","Potencia turboprop [W]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','BG1');

%% Nº rotor, Potencia motor eletrico
clearvars -except constants
data = load_project('transport.json');
X_values = [2 4 6 8];
Y_values = linspace(150000,600000,15);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.vehicle.components{12, 1}.number = X_values(i);
        data.vehicle.components{16, 1}.number = X_values(i);
        data.vehicle.components{16, 1}.max_power = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Nº rotors","Potencia motor eletrico [W]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','BO1');

%% Angulo subida, velocidade subida
clearvars -except constants
data = load_project('transport.json');
X_values = linspace(1,9,15);
Y_values = linspace(40,90,15);
j = 1;

for i=1:length(X_values)
    for ii=1:length(Y_values)
        data.mission.segments{5, 1}.angle = X_values(i);
        data.mission.segments{4, 1}.velocity(2) = Y_values(ii);
        data.mission.segments{5, 1}.velocity = Y_values(ii);
        
        %adt
        data.mission = build_mission(data.mission);
        data.vehicle = build_vehicle(data.mission, data.vehicle);
        data.vehicle = aero_analysis(data.mission, data.vehicle);
        [data.mission, data.vehicle] = mass_analysis(data.mission, data.vehicle, data.energy);
        [WS(j,1), WA(j,1), WPf(j,1), WPv(j,1)] = Auxxx_imp_param(data.mission, data.vehicle, data.energy);
        X(j,1) = X_values(i);
        Y(j,1) = Y_values(ii);
        MTOW(j,1) = data.vehicle.mass;
        j = j + 1;
    end
end

%table
varNames = ["Angulo subida [deg]","Velocidade subida [m/s]","MTOW [kg]","W/S [N/m^2]","W/A [N/m^2]","W/Pf [N/W]","W/Pv [N/W]"];
tab = table(X,Y,MTOW,WS,WA,WPf,WPv,'VariableNames',varNames);
writetable(tab,'VarData.xlsx','Sheet',1,'Range','BW1');