function[tracernew,bottomFlux]=mixing(tracer,z,z_edge,Kz,dt,deep)

% deep=Deep;
% z=z(:,3);
% z_edge=z_edge(:,3);
% Kz=Kz_edge(:,4);
% dt=dt_phys;
% clear Flux numlevs deriv

numlevs=length(tracer(:,1));
deltaz=z(2:end)-z(1:end-1);
deltaz_edge=z_edge(2:end)-z_edge(1:end-1);



tracer(numlevs+1,:)=deep;
tracernew=tracer;

Flux=zeros(size(tracernew(1:end-1,:)));
for i=1:length(tracernew(1,:))
    deriv=(tracer(1:end-1,i)-tracer(2:end,i))./deltaz;
    Flux(:,i)=-deriv.*Kz(2:end-1)*dt;   
end
bottomFlux=Flux(end,:);

for i=1:length(tracernew(1,:))
    tracernew(2:end,i)=tracernew(2:end,i)-Flux(1:end,i)./deltaz_edge(2:end);
    tracernew(1:end-1,i)=tracernew(1:end-1,i)+Flux(1:end,i)./deltaz_edge(1:end-1);
end

tracernew=tracernew(1:end-1,:);

