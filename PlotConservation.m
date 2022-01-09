function [output] = PlotConservation(tracer,tracer_init,z,t,N15on,iN15,Thoriumon,iThorium,euponly)

if euponly==1
    'You will not get conservation in the euphotic zone only version of the model'
end

figure(9)
set(gcf,'Position',[200 100 1400 600])

subplot(1,2+Thoriumon+N15on,1)
plot(sum(tracer(:,1:13)')'./sum(tracer_init(1:end-1,4:16)')',z(1:end-1,3),'-k','LineWidth',2,'Color','b')
ylabel('Depth (m)')
xlabel('Nitrogen (Ratio to initial)')
title(['time = ',num2str(t),' days'])
set(gca,'box','on')
axis ij
set(gcf,'Color','w')

subplot(1,2+Thoriumon+N15on,2)
r_SI_N = 1.0;
plot((r_SI_N*tracer(:,2)'+tracer(:,14)'+tracer(:,15)'+tracer(:,16)')./(r_SI_N*tracer_init(1:end-1,5)'+tracer_init(1:end-1,17)'+tracer_init(1:end-1,18)'+tracer_init(1:end-1,19)'),z(1:end-1,3),'-k','LineWidth',2,'Color','b')
ylabel('Depth (m)')
xlabel('Silica (Ratio to initial)')
title(['time = ',num2str(t),' days'])
set(gca,'box','on')
axis ij
set(gcf,'Color','w')

if N15on==1
    subplot(1,2+Thoriumon+N15on,3)
    plot(sum(tracer(:,iN15:iN15+12)')'./sum(tracer_init(1:end-1,iN15+3:iN15+12+3)')',z(1:end-1,3),'-k','LineWidth',2,'Color',[0.2 1 0.2])
    ylabel('Depth (m)')
    xlabel('^1^5N (Ratio to initial)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
    pause(0.2)
end

if Thoriumon==1
    subplot(1,2+Thoriumon+N15on,3+N15on)
    plot(sum(tracer(:,iThorium:iThorium+9)')'./sum(tracer_init(1:end-1,iThorium+3:iThorium+9+3)')',z(1:end-1,3),'-k','LineWidth',2,'Color',[0.2 1 0.2])
    ylabel('Depth (m)')
    xlabel('Total Thorium (Ratio to initial)')
    title(['time = ',num2str(t),' days'])
    set(gca,'box','on')
    axis ij
    set(gcf,'Color','w')
    pause(0.2)
end

output=0;