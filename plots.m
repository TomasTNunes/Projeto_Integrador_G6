clear all
close all
clc

filename = 'Dados PAer.xlsx';
fig_n = 1;

%% raio rotores, numero rotores, W/A
clearvars -except fig_n filename
xlRange = 'B8:F19';
xlRange_var = 'B7:F7';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
data_Z = data(:,5);
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
zlabel(data_varname{5})
fig_n = fig_n + 1;


%% h vertical climb, velocidade, W/S
clearvars -except fig_n filename
xlRange = 'K8:O31';
xlRange_var = 'K7:O7';
[~,data_varname] = xlsread(filename,xlRange_var);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
data_Z = data(:,5);
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
zlabel(data_varname{5})
fig_n = fig_n + 1;
% plot(data_Y(1:6),data_Z(1:6))
% hold on
% plot(data_Y(7:12),data_Z(7:12))
% plot(data_Y(13:18),data_Z(13:18))
% plot(data_Y(19:24),data_Z(19:24))
