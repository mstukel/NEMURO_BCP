function[Coeff0,Coeff1,Coeff2,BottomCoeff]=CalculateMixingCoefficients(tracer,z,z_edge,Kz,dt,deep)

deltaz=z(2:end)-z(1:end-1);
deltaz_edge=z_edge(2:end)-z_edge(1:end-1);

tracer(end+1,:)=deep;


for j=1:length(tracer(1,:))
    i=1;
    Coeff0(i,j) = 1 - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
    Coeff1(i,j) = 0;
    Coeff2(i,j) = (Kz(i+1)*dt)./(deltaz_edge(i)*deltaz(i));
    for i=2:length(tracer(:,1))-1
        i;
        Coeff0(i,j) = 1 - (Kz(i) * dt)./(deltaz_edge(i).*deltaz(i-1)) - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
        Coeff1(i,j) = (Kz(i)*dt)./(deltaz_edge(i)*deltaz(i-1));
        Coeff2(i,j) = (Kz(i+1)*dt)./(deltaz_edge(i)*deltaz(i));
    end
end


% tracernew = tracer(1:end-1,:).*Coeff0 + ...
%     [zeros(1,length(tracer(1,:)));tracer(1:end-2,:)].*Coeff1 + ...
%     tracer(2:end,:).*Coeff2;

BottomCoeff = Kz(end-1)*dt/deltaz(end);

% bottomflux=(tracer(end,:)-tracer(end-1,:))*BottomCoeff;
