clear all
close all
clc

filename = 'Centros_Geometricos.xlsx';

xlRange = 'B2:D15';
xlRange_name = 'A2:A15';
[~,data_partname] = xlsread(filename,xlRange_name);
data = xlsread(filename,xlRange);
data_X = data(:,1);
data_Y = data(:,2);
data_Z = data(:,3);


data_X_corr = data_X - data_X(1);
data_Y_corr = data_Y - data_Y(1);
data_Z_corr = data_Z - data_Z(1);


R = [-1  0  0;
      0  0 -1;
      0 -1  0;];
for i=1:length(data_X)
    pos(:,i) = R*[data_X_corr(i);data_Y_corr(i);data_Z_corr(i)];
end
X = pos(1,:);
Y = pos(2,:);
Z = pos(3,:);

plotColors = jet(length(X));

figure()
for ii=1:length(X)
    if ii==1
        plot3(X(ii),Y(ii),Z(ii),'xk','LineWidth', 1.5)
    else
        plot3(X(ii),Y(ii),Z(ii),'o','Color', plotColors(ii,:), 'LineWidth', 1.5)
    end
    grid on
    hold on
end

xlabel('X')
ylabel('Y')
zlabel('Z')
set ( gca, 'ydir', 'reverse' )
set ( gca, 'zdir', 'reverse' )

legend(data_partname)