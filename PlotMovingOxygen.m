function [output] = PlotMovingOxygen(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15)

numpanels = 3;

filename=['zz.Cycle',num2str(Cycle(1)),'-',num2str(Cycle(2)),' Oxygen.gif'];

    figure(538)
    clf
    set(gcf,'Position',[100 50 1400 700])
    subplot(1,numpanels,1)
    hold on
    %plot(tracer(:,19),z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
    %plot(tracer(:,20),z(1:end-1,3),'-sk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.4 0.1])
    plot(tracer(:,19),z(1:end-1,3),'-k','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0.5])
    %legend('DIC','Alk','Oxy','Location','SouthEast')
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('Oxygen (mmol O m^-^3)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
    xlim([0 350])

if CO2on==1

    
    
    
    subplot(1,numpanels,2)
    hold on
    plot(tracer(:,20),z(1:end-1,3),'-k','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
    plot(tracer(:,21),z(1:end-1,3),'-k','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.4 0.1])
    %plot(tracer(:,21),z(1:end-1,3),'-^k','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0.5])
    legend('DIC','Alk','Location','SouthEast')
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('Carbon System (mmol C m^-^3)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
    xlim([0 4000])
    
    
    
    subplot(1,numpanels,3)
    hold on
    plot(sum(tracer(:,[1:7,10:13])')',z(1:end-1,3),'-k','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
    plot(sum(tracer(:,[8:9])')',z(1:end-1,3),'-k','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
    %plot(tracer(:,21),z(1:end-1,3),'-^k','Color',[0.5 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0.5])
    legend('Org N','Nuts','Location','SouthEast')
    set(gca,'FontSize',12)
    ylabel('Depth (m)')
    xlabel('Nitrogen (mmol N m^-^3)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
    xlim([0 50])
end




pause(0.1)
% frame = getframe(2);
% mov = frame2im(frame);
% [imind,cm] = rgb2ind(mov,256);
% if round(t)==1
%     imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.2,'TransparentColor',-1);
% else
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2,'TransparentColor',-1);
% end

output=0;