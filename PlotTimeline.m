function [output] = PlotTimeline(t,tracer,z)

figure(1)
set(gcf,'Position',[50 50 1200 600])
subplot(1,2,1)
plot(t,tracer(1,1),'.g')
hold on
plot(t,tracer(1,2),'.g','MarkerEdgeColor',[0 0.5 0])
plot(t,tracer(1,3),'.b')
plot(t,tracer(1,4),'.r','MarkerEdgeColor',[1 0 0])
plot(t,tracer(1,5),'.r','MarkerEdgeColor',[1 0 1])
plot(t,tracer(1,6),'.r','MarkerEdgeColor',[0.5 0 0])
plot(t,tracer(1,7),'.r','MarkerEdgeColor',[0.5 0 0.5])
title(['Depth = ',num2str(z(1,3))])
plot(t,sum(tracer(1,1:13)),'.k')
%legend('SP','LP','SZ','LZres','LZdvm','PZres','PZdvm')


subplot(1,2,2)
plot(t,tracer(8,1),'.g')
hold on
plot(t,tracer(8,2),'.g','MarkerEdgeColor',[0 0.5 0])
plot(t,tracer(8,3),'.b')
plot(t,tracer(8,4),'.r','MarkerEdgeColor',[1 0 0])
plot(t,tracer(8,5),'.r','MarkerEdgeColor',[1 0 1])
plot(t,tracer(8,6),'.r','MarkerEdgeColor',[0.5 0 0])
plot(t,tracer(8,7),'.r','MarkerEdgeColor',[0.5 0 0.5])
title(['Depth = ',num2str(z(8,3))])
plot(t,sum(tracer(8,1:13)),'.k')
%legend('SP','LP','SZ','LZres','LZdvm','PZres','PZdvm')

output=0;