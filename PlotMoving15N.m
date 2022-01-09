function [output]=PlotMoving15N(Cycle,tracer,z,t,iN15)


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
RPON = tracer(:,iN15+9)./tracer(:,10);
RLPON = tracer(:,iN15+10)./tracer(:,11);
dsp = (Rsp-RN2)/RN2*1000;
dlp = (Rlp-RN2)/RN2*1000;
dsz = (Rsz-RN2)/RN2*1000;
dlzres = (Rlzres-RN2)/RN2*1000;
dlzdvm = (Rlzdvm-RN2)/RN2*1000;
dpzres = (Rpzres-RN2)/RN2*1000;
dpzdvm = (Rpzdvm-RN2)/RN2*1000;
dNO3 = (RNO3-RN2)/RN2*1000;
dNH4 = (RNH4-RN2)/RN2*1000;
dPON = (RPON-RN2)/RN2*1000;
dLPON = (RLPON-RN2)/RN2*1000;

figure(94)
set(gcf,'Position',[50 50 1000 600])
subplot(1,4,1)
hold off
plot(tracer(:,5),z(1:end-1,3),'-r')
hold on
plot(tracer(:,7),z(1:end-1,3),'-.r','Color',[0.5 0 0])
plot(tracer(:,8),z(1:end-1,3),'-.b','Color',[0 0 0.5])
plot(tracer(:,9),z(1:end-1,3),'-b','Color',[0.2 0.2 1])
plot(tracer(:,10),z(1:end-1,3),'-m','Color',[1 0.8 0.05])
plot(tracer(:,11),z(1:end-1,3),'-m','Color',[0.7 0.5 0])
xlabel('LZdvm and PZdvm')
legend('LZdvm','PZdvm','NO3','NH4','PON','LPON','Location','SouthEast')
axis ij

subplot(1,4,2)
hold off
plot(tracer(:,iN15+4),z(1:end-1,3),'-r')
hold on
plot(tracer(:,iN15+6),z(1:end-1,3),'-.r','Color',[0.5 0 0])
plot(tracer(:,iN15+7),z(1:end-1,3),'-.b','Color',[0 0 0.5])
plot(tracer(:,iN15+8),z(1:end-1,3),'-b','Color',[0.2 0.2 1])
plot(tracer(:,iN15+9),z(1:end-1,3),'-m','Color',[1 0.8 0.05])
plot(tracer(:,iN15+10),z(1:end-1,3),'-m','Color',[0.7 0.5 0])
xlabel('^1^5N (absolute value)')
axis ij

subplot(1,4,3)
hold off
plot(Rlzdvm,z(1:end-1,3),'-r')
hold on
plot(Rpzdvm,z(1:end-1,3),'-.r','Color',[0.5 0 0])
plot(RNO3,z(1:end-1,3),'-.b','Color',[0 0 0.5])
plot(RNH4,z(1:end-1,3),'-b','Color',[0.2 0.2 1])
plot(RPON,z(1:end-1,3),'-m','Color',[1 0.8 0.05])
plot(RLPON,z(1:end-1,3),'-m','Color',[0.7 0.5 0])
xlabel('^1^5N:1^4N (ratio)')
axis ij

subplot(1,4,4)
hold off
plot(dlzdvm,z(1:end-1,3),'-r')
hold on
plot(dpzdvm,z(1:end-1,3),'-.r','Color',[0.5 0 0])
plot(dNO3,z(1:end-1,3),'-.b','Color',[0 0 0.5])
plot(dNH4,z(1:end-1,3),'-b','Color',[0.2 0.2 1])
plot(dPON,z(1:end-1,3),'-m','Color',[1 0.8 0.05])
plot(dLPON,z(1:end-1,3),'-m','Color',[0.7 0.5 0])
xlabel('d^1^5N')
axis ij
pause(0.5)
