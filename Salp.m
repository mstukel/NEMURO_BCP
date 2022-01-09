function [tracer_new,export,exportsi] = Salp(tracer,clearancerate,exportfrac,dt)

frac2nh4 = 0.4;
frac2donlab = 0.3;
frac2ponsmall = 0.3;

tracer_new = tracer;

tracer_new(:,[1:3,17:18])=tracer(:,[1:3,17:18])*(1-clearancerate*dt);  %salps feed on PS, PL, and ZS, and also remove Chl
consumed=tracer(:,1)*clearancerate*dt+tracer(:,2)*clearancerate*dt+tracer(:,3)*clearancerate*dt;
consumedsi = tracer(:,2)*clearancerate*dt;
for i=1:length(tracer(:,1))
    export(i,1) = sum(consumed(1:i))*exportfrac;
    exportsi(i,1) = sum(consumedsi(1:i))*exportfrac;
end

tracer_new(:,9)=tracer(:,9)+consumed*(1-exportfrac)*frac2nh4;
tracer_new(:,10)=tracer(:,10)+consumed*(1-exportfrac)*frac2donlab;
tracer_new(:,12)=tracer(:,12)+consumed*(1-exportfrac)*frac2ponsmall;

tracer_new(:,14)=tracer(:,14)+consumedsi*(1-exportfrac)*(frac2nh4+frac2donlab);
tracer_new(:,15)=tracer(:,15)+consumedsi*(1-exportfrac)*(frac2ponsmall);
