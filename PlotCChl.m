function [output] = PlotCChl(Cycle,tracer,PAR,z,t,dt)


filename=['zz.Cycle',num2str(Cycle(1)),'-',num2str(Cycle(2)),' CChl.gif'];
figure(3)
set(gcf,'Position',[50 50 1400 600])
clf
subplot(1,4,1)
PS_mgC = tracer(:,1)*6.625*12;
PL_mgC = tracer(:,2)*6.625*12;
plot(PS_mgC,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(PL_mgC,z(1:end-1,3),'-ok','Color',[0 0 0.5],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.5])
legend('PS','PL','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Phytoplankton Biomass (mg C m^-^3)')
axis ij

subplot(1,4,2)
PS_mgC = tracer(:,1)*6.625*12;
PL_mgC = tracer(:,2)*6.625*12;
CChl_PS = PS_mgC./tracer(:,17);
CChl_PL = PL_mgC./tracer(:,18);
plot(CChl_PS,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(CChl_PL,z(1:end-1,3),'-ok','Color',[0 0 0.5],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.5])
legend('PS','PL','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Phytoplankton C:Chl Ratio (g:g)')
axis ij

subplot(1,4,3)
PS_mgC = tracer(:,1)*6.625*12;
PL_mgC = tracer(:,2)*6.625*12;
CChl_PS = PS_mgC./tracer(:,17);
CChl_PL = PL_mgC./tracer(:,18);
plot(1./CChl_PS,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(1./CChl_PL,z(1:end-1,3),'-ok','Color',[0 0 0.5],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0.5])
legend('PS','PL','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Phytoplankton Chl:C (g:g)')
axis ij

subplot(1,4,4)
plot(log10(PAR),z(1:end-1,3),'-k','Color',[0.8 0.8 0.2],'LineWidth',2)
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('PAR')
set(gca,'XTick',log10([10^-3:10^-3:10^-2,2*10^-2:10^-2:0.1,0.2:0.1:1,2:10,20:10:100,200:100:1000]))
set(gca,'XTickLabel',{'10^-^3','','','','','','','','','10^-^2','','','','','','','','',...
    '10^-^1','','','','','','','','','1','','','','','','','','',...
    '10','','','','','','','','','100','','','','','','','','','1000'})
axis ij

pause(0.1)

frame = getframe(3);
mov = frame2im(frame);
[imind,cm] = rgb2ind(mov,256);
if round(t)==1
    imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.2,'TransparentColor',-1);
else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2,'TransparentColor',-1);
end


output=0;