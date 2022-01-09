function[Coeff0,Coeff1,Coeff2,swimspd]=CalculateDVMCoefficients(tracer,z,z_edge,dt,targetdepth)

%--------------------------------------------------
%Hard-wired diel vertical migration swimming speeds
%--------------------------------------------------
Kz = 10^-1 * 24 * 60 * 60;
maxswim = 0.03 * 60 * 60 * 24;    %30 cm / s  (Bianchi et al. 2013)
minswim = 1 * 24;     %2 m/h 
tuningparam = 200;      %This parameter determines how quickly the zooplankton swimming speed transitions from max swimspeed to minswimspd
%--------------------------------------------------
%--------------------------------------------------


%--------------------------------------------------
%Preparatory calculations--------------------------
%--------------------------------------------------
if targetdepth==round(targetdepth)
    targetdepth=targetdepth+0.1;   %This code is simply to avoid divide by zeros that can occur if the targetdepth is exactly equal to the midpoint of a depth bin
end
swimspd = (z - targetdepth)./abs(z - targetdepth) .* ...
    (minswim + (maxswim-minswim) .* (abs(z-targetdepth)./((abs(z-targetdepth)+tuningparam))));
deltaz=z(2:end)-z(1:end-1);
deltaz_edge=z_edge(2:end)-z_edge(1:end-1);
%--------------------------------------------------
%--------------------------------------------------


i=1;
Coeff0(i,1) = 1 - (Kz *dt)./(deltaz_edge(i)*deltaz(i)) - max(0,-swimspd(i)*dt./deltaz_edge(i));
Coeff1(i,1) = 0;
Coeff2(i,1) = (Kz*dt)./(deltaz_edge(i)*deltaz(i)) + max(0,swimspd(i+1)*dt./deltaz_edge(i));
for i=2:length(tracer(:,1))-1
    Coeff0(i,1) = 1 - (Kz * dt)./(deltaz_edge(i).*deltaz(i-1)) - (Kz *dt)./(deltaz_edge(i)*deltaz(i)) - ...
        max(0,swimspd(i)*dt./deltaz_edge(i)) - max(0,-swimspd(i)*dt./deltaz_edge(i));
    Coeff1(i,1) = (Kz*dt)./(deltaz_edge(i)*deltaz(i-1)) + max(0,-swimspd(i-1)*dt./deltaz_edge(i));
    Coeff2(i,1) = (Kz*dt)./(deltaz_edge(i)*deltaz(i)) + max(0,swimspd(i+1)*dt./deltaz_edge(i));
end
i=length(tracer(:,1));
Coeff0(i,1) = 1 - (Kz * dt)./(deltaz_edge(i).*deltaz(i-1)) - (max(0,swimspd(i))*dt)./deltaz_edge(i-1);
Coeff1(i,1) = (Kz*dt)./(deltaz_edge(i)*deltaz(i-1)) + max(0,-swimspd(i-1)*dt./deltaz_edge(i-1));
Coeff2(i,1) = 0;

for j=2:length(tracer(1,:))
    Coeff0(:,j)=Coeff0(:,1);
    Coeff1(:,j)=Coeff1(:,1);
    Coeff2(:,j)=Coeff2(:,1);
end


% tracernew = tracer(1:end-1,:).*Coeff0 + ...
%     [zeros(1,length(tracer(1,:)));tracer(1:end-2,:)].*Coeff1 + ...
%     tracer(2:end,:).*Coeff2;
 