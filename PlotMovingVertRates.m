function [output]=PlotMovingVertRates(Cycle,tracer,GPP,NPP,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,mu_sp,mu_lp,m_chl,mu_chl,NO3up,NH4up,Siup,z,t,dt)


filename=['zz.Cycle',num2str(Cycle(1)),'-',num2str(Cycle(2)),' Rates.gif'];
figure(3)
set(gcf,'Position',[50 50 1400 600])
clf
subplot(1,3,1)
plot(GPP*6.625*12/dt,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(NPP*6.625*12/dt,z(1:end-1,3),'-ok','Color',[0 0.5 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot((SP2SZ+LP2SZ)*6.625*12/dt,z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1])
plot((LZresgraz+LZdvmgraz)*6.625*12/dt,z(1:end-1,3),'-sk','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
plot((PZresgraz+PZdvmgraz)*6.625*12/dt,z(1:end-1,3),'-sk','Color',[1 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1])
legend('GPP','NPP','Protist Grazing','LZ Grazing','PZ Grazing','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Rate (mg C m^-^3 d^-^1)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,3,2)
plot(mu_sp/dt,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(mu_lp/dt,z(1:end-1,3),'-ok','Color',[0 0.5 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot(m_chl/dt,z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
plot(mu_chl/dt,z(1:end-1,3),'-ok','Color',[1 0 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.4 0.1])
plot(NO3up./(NO3up+NH4up),z(1:end-1,3),'--k','Color',[0.5 0.5 0.5],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
legend('mu_S_P','mu_L_P','m_c_h_l','mu_c_h_l','f-ratio','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Rate (d^-^1) or Ratio (unitless)')
title(['time = ',num2str(t),' days'])
axis ij

subplot(1,3,3)
plot(NH4up,z(1:end-1,3),'-ok','Color','g','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(NO3up,z(1:end-1,3),'-ok','Color',[0 0.5 0],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0])
plot(Siup,z(1:end-1,3),'-dk','Color',[0 0 1],'LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
legend('NH4_u_p','NO3_u_p','Si_u_p','Location','SouthEast')
set(gca,'FontSize',12)
ylabel('Depth (m)')
xlabel('Rate (mmol m^-^3 d^-^1)')
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