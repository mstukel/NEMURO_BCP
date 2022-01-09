function [output] = PlotMovingVertical(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15)

numpanels = 3 + CO2on + Thoriumon + N15on;

filename=['zz.Cycle',num2str(Cycle(1)),'-',num2str(Cycle(2)),' Tracers.gif'];
figure(2)
set(gcf,'Position',[50 50 1800 600])
set(gca,'Color','w')
clf
subplot(1,numpanels,1)
plot(tracer(:,1),z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(tracer(:,2),z(1:end-1,3),'-ok','Color',[0 0.5 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot(tracer(:,3),z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
plot(tracer(:,4),z(1:end-1,3),'-sk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
plot(tracer(:,5),z(1:end-1,3),'-sk','Color',[1 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1])
plot(tracer(:,6),z(1:end-1,3),'-sk','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0 0])
plot(tracer(:,7),z(1:end-1,3),'-sk','Color',[0.5 0 0.5],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0 0.5])
legend('SP','LP','SZ','LZres','LZdvm','PZres','PZdvm','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Living Tracers (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,numpanels,2)
plot(tracer(:,8),z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(tracer(:,9),z(1:end-1,3),'-ok','Color',[0 0.5 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot(tracer(:,14),z(1:end-1,3),'-^k','Color',[1 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
legend('NO3','NH4','Si','Location','NorthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Nutrient Tracers (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,numpanels,3)
hold on
plot(tracer(:,10),z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
plot(tracer(:,12),z(1:end-1,3),'-sk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.4 0.1])
plot(tracer(:,15),z(1:end-1,3),'-^k','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0.5])
legend('PON','DON_l_a_b','Opal','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Non-living Tracers (mmol N m^-^3)')
title(['time = ',num2str(t),' days'])
set(gca,'box','on')
axis ij

if CO2on==1
    subplot(1,numpanels,4)
    hold on
    plot(tracer(:,20),z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
    plot(tracer(:,21),z(1:end-1,3),'-sk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.4 0.1])
    plot(tracer(:,19),z(1:end-1,3),'-^k','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0.5])
    legend('DIC','Alk','Oxy','Location','SouthEast')
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('Carbon system (mmol C m^-^3)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
end

if Thoriumon==1
    subplot(1,numpanels,4+CO2on)
    U238 = salinity(1:length(tracer(:,1))).*0.0786-0.315;
    hold on
    plot(U238,z(1:end-1,3),'-k','LineWidth',3,'Color',[0.2 0.2 0.2])
    plot(tracer(:,22),z(1:end-1,3),'-dk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
    plot(sum(tracer(:,23:31)')',z(1:end-1,3),'-sk','Color',[0 1 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 1 0])
    plot(sum(tracer(:,22:31)')',z(1:end-1,3),'-k','LineWidth',3,'Color',[0 0 0.8])
    legend('^2^3^8U','Dissolved ^2^3^4Th','Particulate ^2^3^4Th','Total ^2^3^4Th','Location','SouthEast')
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('Activity (dpm L^-^1)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
end

if N15on==1
    subplot(1,numpanels,4+CO2on+Thoriumon)
    RN2=0.0036765;  
    Rsp = tracer(:,iN15)./tracer(:,1);
    Rlp = tracer(:,iN15+1)./tracer(:,2);
    Rsz = tracer(:,iN15+2)./tracer(:,3);
    Rlzres = tracer(:,iN15+3)./tracer(:,4);
    Rlzdvm = tracer(:,iN15+4)./tracer(:,5);
    Rpzres = tracer(:,iN15+5)./tracer(:,6);
    Rpzdvm = tracer(:,iN15+6)./tracer(:,7);
    RNO3 = tracer(:,iN15+7)./tracer(:,8);
    RNH4 = tracer(:,iN15+8)./tracer(:,9);
    dsp = (Rsp-RN2)/RN2*1000;
    dlp = (Rlp-RN2)/RN2*1000;
    dsz = (Rsz-RN2)/RN2*1000;
    dlzres = (Rlzres-RN2)/RN2*1000;
    dlzdvm = (Rlzdvm-RN2)/RN2*1000;
    dpzres = (Rpzres-RN2)/RN2*1000;
    dpzdvm = (Rpzdvm-RN2)/RN2*1000;
    dNO3 = (RNO3-RN2)/RN2*1000;
    dNH4 = (RNH4-RN2)/RN2*1000;
    
    hold on
    plot(dNO3,z(1:end-1,3),'-k','LineWidth',2,'Color',[0 0 0.8])
    plot(dNH4,z(1:end-1,3),'-k','LineWidth',2,'Color',[0.3 0.5 1])
    plot(dsp,z(1:end-1,3),'-k','LineWidth',2,'Color',[0.2 1 0.2])
    plot(dlp,z(1:end-1,3),'-k','LineWidth',2,'Color',[0 0.7 0])
    plot(dsz,z(1:end-1,3),'-k','LineWidth',2,'Color','m')
    plot(dlzres,z(1:end-1,3),'-k','LineWidth',2,'Color','r')
    plot(dlzdvm,z(1:end-1,3),'-.k','LineWidth',2,'Color','r')
    plot(dpzres,z(1:end-1,3),'-k','LineWidth',2,'Color',[0.6 0 0])
    plot(dpzdvm,z(1:end-1,3),'-.k','LineWidth',2,'Color',[0.6 0 0])
    h=legend('NO3','NH4','SP','LP','SZ','LZres','LZdvm','PZres','PZdvm','Location','SouthWest');
    set(h,'FontSize',8)
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('d^1^5N')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
end
    



pause(0.1)
frame = getframe(2);
mov = frame2im(frame);
[imind,cm] = rgb2ind(mov,256);
if round(t)==1
    imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.2,'TransparentColor',-1);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2,'TransparentColor',-1);
end

output=0;