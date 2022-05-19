clear all
close all
clc

filename = 'VarData.xlsx';
fig_n = 1;

%% Raio Rotores, Nº Rotores
%MTOW
clearvars -except fig_n filename
xlRange = 'A2:G41';
xlRange_var = 'A1:G61';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(Y_values)
    k = find(data_Y==Y_values(i));
    plot(data_X(k),data_Z(k))
    hold on
end
xlabel(data_varname{1})
ylabel(data_varname{Z_ind})
legendStrings = "Nº rotores = " + string(cast(Y_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%W/A
clearvars -except fig_n filename data_varname data
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 5;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(Y_values)
    k = find(data_Y==Y_values(i));
    plot(data_X(k),data_Z(k))
    hold on
end
xlabel(data_varname{1})
ylabel(data_varname{Z_ind})
legendStrings = "Nº rotores = " + string(cast(Y_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%% Aspect Ratio, Mean Chord
%W/S
clearvars -except fig_n filename
xlRange = 'T2:Z226';
xlRange_var = 'T1:Z1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 4;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
for i=1:length(X_values)
    j = 1;
    for ii=1:data_length
        if data_X(ii) == X_values(i)
            X(j,i) = data_X(ii);
            Y(j,i) = data_Y(ii);
            Z(j,i) = data_Z(ii);
            j = j + 1;
        end
    end
end

figure(fig_n)
surf(X,Y,Z)
colorbar
xlabel(data_varname{1})
ylabel(data_varname{2})
zlabel(data_varname{Z_ind})
fig_n = fig_n + 1;

%% Crew pax payload, massa fuselagem
%MTOW
clearvars -except fig_n filename
xlRange = 'AB2:AH401';
xlRange_var = 'AB1:AH1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
for i=1:length(X_values)
    j = 1;
    for ii=1:data_length
        if data_X(ii) == X_values(i)
            X(j,i) = data_X(ii);
            Y(j,i) = data_Y(ii);
            Z(j,i) = data_Z(ii);
            j = j + 1;
        end
    end
end

figure(fig_n)
surf(X,Y,Z)
colorbar
xlabel(data_varname{1})
ylabel(data_varname{2})
zlabel(data_varname{Z_ind})
fig_n = fig_n + 1;

%% massa fuselagem, range
%MTOW
clearvars -except fig_n filename
xlRange = 'AJ2:AP401';
xlRange_var = 'AJ1:AP1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
for i=1:length(X_values)
    j = 1;
    for ii=1:data_length
        if data_X(ii) == X_values(i)
            X(j,i) = data_X(ii);
            Y(j,i) = data_Y(ii);
            Z(j,i) = data_Z(ii);
            j = j + 1;
        end
    end
end

figure(fig_n)
surf(X,Y,Z)
colorbar
xlabel(data_varname{1})
ylabel(data_varname{2})
zlabel(data_varname{Z_ind})
fig_n = fig_n + 1;

%% h vertical climb, Velocidade de subida vertical
%MTOW
clearvars -except fig_n filename
xlRange = 'I2:O226';
xlRange_var = 'I1:O1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
for i=1:length(X_values)
    j = 1;
    for ii=1:data_length
        if data_X(ii) == X_values(i)
            X(j,i) = data_X(ii);
            Y(j,i) = data_Y(ii);
            Z(j,i) = data_Z(ii);
            j = j + 1;
        end
    end
end

figure(fig_n)
surf(X,Y,Z)
colorbar
xlabel(data_varname{1})
ylabel(data_varname{2})
zlabel(data_varname{Z_ind})
fig_n = fig_n + 1;

%% bSFC cruise
%MTOW
clearvars -except fig_n filename
xlRange = 'AR2:AW101';
xlRange_var = 'AR1:AW1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);

figure(fig_n)
plot(data_X,data_Y)
xlabel(data_varname{1})
ylabel(data_varname{2})
title('Velocidade cruise = 80 m/s')
fig_n = fig_n + 1;

%% range cruise, v=cte
%MTOW
clearvars -except fig_n filename
xlRange = 'AY2:BE101';
xlRange_var = 'AY1:BE1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,2);
data_Y = data(:,3);

figure(fig_n)
plot(data_X,data_Y)
xlabel(data_varname{2})
ylabel(data_varname{3})
title('Velocidade cruise = 80 m/s')
fig_n = fig_n + 1;

%% Nº propellers, Potencia turboprop
%MTOW
clearvars -except fig_n filename
xlRange = 'BG2:BM46';
xlRange_var = 'BG1:BM1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(X_values)
    k = find(data_X==X_values(i));
    plot(data_Y(k),data_Z(k))
    hold on
end
xlabel(data_varname{2})
ylabel(data_varname{Z_ind})
legendStrings = "Nº propellers = " + string(cast(X_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%W/Pf
clearvars -except fig_n filename data_varname data
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 6;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(X_values)
    k = find(data_X==X_values(i));
    plot(data_Y(k),data_Z(k))
    hold on
end
xlabel(data_varname{2})
ylabel(data_varname{Z_ind})
legendStrings = "Nº propellers = " + string(cast(X_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%% Nº rotor, Potencia motor eletrico
%MTOW
clearvars -except fig_n filename
xlRange = 'BO2:BU61';
xlRange_var = 'BO1:BU1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(X_values)
    k = find(data_X==X_values(i));
    plot(data_Y(k),data_Z(k))
    hold on
end
xlabel(data_varname{2})
ylabel(data_varname{Z_ind})
legendStrings = "Nº Rotores = " + string(cast(X_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%W/A
clearvars -except fig_n filename data_varname data
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 5;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(X_values)
    k = find(data_X==X_values(i));
    plot(data_Y(k),data_Z(k))
    hold on
end
xlabel(data_varname{2})
ylabel(data_varname{Z_ind})
legendStrings = "Nº Rotores = " + string(cast(X_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%W/Pv
clearvars -except fig_n filename data_varname data
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 7;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
figure(fig_n)
for i=1:length(X_values)
    k = find(data_X==X_values(i));
    plot(data_Y(k),data_Z(k))
    hold on
end
xlabel(data_varname{2})
ylabel(data_varname{Z_ind})
legendStrings = "Nº Rotores = " + string(cast(X_values,'int32'));
legend(legendStrings,'Location','Northeast');
fig_n = fig_n + 1;

%% Densidade da bateria, Massa bateria
clearvars -except fig_n filename
xlRange = 'Q2:R101';
xlRange_var = 'Q1:R1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);

figure(fig_n)
subplot(1,2,1)
plot(data_X,data_Y)
xlabel(data_varname{1})
ylabel(data_varname{2})
subplot(1,2,2)
plot(data_X,data_Y)
xlabel(data_varname{1})
ylabel(data_varname{2})
xlim([200000 2000000])
fig_n = fig_n + 1;

%% Angulo subida, velocidade subida
%MTOW
clearvars -except fig_n filename
xlRange = 'BW2:CC226';
xlRange_var = 'BW1:CC1';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
Z_ind = 3;
data_Z = data(:,Z_ind);
data_length = length(data_X);
X_values = unique(data_X);
Y_values = unique(data_Y);
for i=1:length(X_values)
    j = 1;
    for ii=1:data_length
        if data_X(ii) == X_values(i)
            X(j,i) = data_X(ii);
            Y(j,i) = data_Y(ii);
            Z(j,i) = data_Z(ii);
            j = j + 1;
        end
    end
end

figure(fig_n)
surf(X,Y,Z)
colorbar
xlabel(data_varname{1})
ylabel(data_varname{2})
zlabel(data_varname{Z_ind})
fig_n = fig_n + 1;
