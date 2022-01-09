function [output] = PlotMovingDVM(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15,dt)

numpanels = 3;

filename=['zz.Cycle',num2str(Cycle(1)),'-',num2str(Cycle(2)),' DVM.gif'];
figure(2)
set(gcf,'Position',[50 50 1800 600])
set(gca,'Color','w')
clf
subplot(1,numpanels,1)
plot(tracer(:,5),z(1:end-1,3),'-ok','Color','r','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
hold on
plot(tracer(:,7),z(1:end-1,3),'-ok','Color',[0.5 0 0],'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0 0])
legend('dvm LZ','dvm PZ','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Vertical Migrators (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,numpanels,2)
plot(tracer(:,4),z(1:end-1,3),'-ok','Color','r','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r')
hold on
plot(tracer(:,6),z(1:end-1,3),'-ok','Color',[0.5 0 0],'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0 0])
legend('res LZ','res PZ','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Non-vertical migrators (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,numpanels,3)
hold on
plot(tracer(:,1),z(1:end-1,3),'-dk','Color',[0 1 0],'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
plot(tracer(:,2),z(1:end-1,3),'-sk','Color',[0 0.5 0],'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot(tracer(:,3),z(1:end-1,3),'-^k','Color',[0.5 0 0.5],'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0 0.5])
legend('SP','LP','SZ','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Phytoplankton (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
set(gca,'box','on')
axis ij




pause(0.1)
frame = getframe(2);
mov = frame2im(frame);
[imind,cm] = rgb2ind(mov,256);
if t==dt
    imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.2,'TransparentColor',-1);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2,'TransparentColor',-1);
end

output=0;