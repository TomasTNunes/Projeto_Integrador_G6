function plot_geo_centers(X,Y,Z,comp)

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

xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
set ( gca,'ydir','reverse')
set ( gca,'zdir','reverse')
title('Centros Geometricos das Componentes')
legend(comp)

end