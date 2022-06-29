function plot_centers(X,Y,Z,comp)

w_CG = 0;
plotColors = jet(length(X));
figure()
for ii=1:length(X)
    if ii==1
        plot3(X(ii),Y(ii),Z(ii),'xk','LineWidth', 1.5)
    elseif strcmp(comp(ii),'CG')
        plot3(X(ii),Y(ii),Z(ii),'*k','LineWidth', 1.8)
        w_CG = 1;
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
if w_CG
    title('Centros Geometricos das Componentes e CG')
else
    title('Centros Geometricos das Componentes')
end
legend(comp)

end