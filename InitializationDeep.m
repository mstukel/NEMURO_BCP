function [tracer,tracer_init,temperature,salinity,Deep,tmpdeep,Kz_mid,Kz_edge,z,z_edge,z_thick,PAR_surf,wind,iCO2,iThorium,iN15,fastsinkingIndices,slowsinkingIndices,mixedIndices,lzdvmIndices,pzdvmIndices,tot_iter]=InitializationDeep(Cycle,CO2on,Thoriumon,N15on,sim_dur,dt)


tot_iter=sim_dur/dt;       % total iterations 

iCO2 = 20;                           %Index of first CO2 submodule variable
iThorium = 20+2*CO2on;               %Index of first Thorium submodule variable
iN15 = 20+2*CO2on+10*Thoriumon;      %Index of first N15 submodule variable



%Loading in initial conditions and forcing
load('Initial Conditions Deep.mat')
% load('..\..\Initial Conditions\Initial Carbon-Nutrients from GLODAP.mat')
% DIC=DIC(find(DIC(:,1)==Cycle(1)+Cycle(2)/1000),:);
% Alk=Alk(find(Alk(:,1)==Cycle(1)+Cycle(2)/1000),:);
% oxy=oxy(find(oxy(:,1)==Cycle(1)+Cycle(2)/1000),:);
% GLODAP_sal=GLODAP_sal(find(GLODAP_sal(:,1)==Cycle(1)+Cycle(2)/1000),:);
load('VerticalDiffusivityAtGridPoints.mat')
load('Temperature-Salinity Deep.mat')
tracer_init=tracer_init(find(tracer_init(:,1)==Cycle(1) & abs(tracer_init(:,2)-Cycle(2))<10^-6),:);
salinity=salinity(find(salinity(:,1)==Cycle(1) & salinity(:,2)==Cycle(2)),4);
temperature=temperature(find(temperature(:,1)==Cycle(1) & temperature(:,2)==Cycle(2)),:);


% for i=1:length(tracer_init(:,1))
%     ind=find(oxy(:,2)==tracer_init(i,3));
%     tracer_init(i,iCO2+3+0)=DIC(ind,3);    %Note that the +3 is because tracer_init has 3 columns with metadata
%     tracer_init(i,iCO2+3+1)=Alk(ind,3);
%     tracer_init(i,iCO2+3+2)=oxy(ind,3);
%     salinity(i,1)=GLODAP_sal(ind,3);
% end
if Thoriumon==1
    CTh = 10; %umol C / dpm;
    U238=salinity*0.0786-0.315;
    tracer_init(:,iThorium+3+1)=tracer_init(:,1+3)*6.625/CTh;
    tracer_init(:,iThorium+3+2)=tracer_init(:,2+3)*6.625/CTh;
    tracer_init(:,iThorium+3+3)=tracer_init(:,3+3)*6.625/CTh;
    tracer_init(:,iThorium+3+4)=tracer_init(:,4+3)*6.625/CTh/10;
    tracer_init(:,iThorium+3+5)=tracer_init(:,5+3)*6.625/CTh/10;
    tracer_init(:,iThorium+3+6)=tracer_init(:,6+3)*6.625/CTh/10;
    tracer_init(:,iThorium+3+7)=tracer_init(:,7+3)*6.625/CTh/10;
    tracer_init(:,iThorium+3+8)=tracer_init(:,10+3)*6.625/CTh;
    tracer_init(:,iThorium+3+9)=tracer_init(:,11+3)*6.625/CTh;
    tracer_init(:,iThorium+3+0)=max(U238-sum(tracer_init(:,23+3:31+3)')',0.1);
end
if N15on==1
    del15NdeepNO3=4;
    RN2=0.0036765;
    RDeepNO3 = del15NdeepNO3/1000*RN2+RN2;
    tracer_init(:,iN15+3+0)=tracer_init(:,3+1)*RDeepNO3;
    tracer_init(:,iN15+3+1)=tracer_init(:,3+2)*RDeepNO3;
    tracer_init(:,iN15+3+2)=tracer_init(:,3+3)*RDeepNO3;
    tracer_init(:,iN15+3+3)=tracer_init(:,3+4)*RDeepNO3;
    tracer_init(:,iN15+3+4)=tracer_init(:,3+5)*RDeepNO3;
    tracer_init(:,iN15+3+5)=tracer_init(:,3+6)*RDeepNO3;
    tracer_init(:,iN15+3+6)=tracer_init(:,3+7)*RDeepNO3;
    tracer_init(:,iN15+3+7)=tracer_init(:,3+8)*RDeepNO3;
    tracer_init(:,iN15+3+8)=tracer_init(:,3+9)*RDeepNO3;
    tracer_init(:,iN15+3+9)=tracer_init(:,3+10)*RDeepNO3;
    tracer_init(:,iN15+3+10)=tracer_init(:,3+11)*RDeepNO3;
    tracer_init(:,iN15+3+11)=tracer_init(:,3+12)*RDeepNO3;
    tracer_init(:,iN15+3+12)=tracer_init(:,3+13)*RDeepNO3;
end

    
tracer=tracer_init(:,4:end);
tracer(:,17:18)=tracer(:,1:2);
Kz_mid=Kz_mid(find(Kz_mid(:,1)==Cycle(1) & Kz_mid(:,2)==Cycle(2)),:);
Kz_mid(:,4)=Kz_mid(:,4)*86400;
Kz_edge=Kz_edge(find(Kz_edge(:,1)==Cycle(1) & Kz_edge(:,2)==Cycle(2)),:);
Kz_edge(:,4)=Kz_edge(:,4)*86400;
z=z(find(z(:,1)==Cycle(1) & abs(z(:,2)-Cycle(2))<10^-6),:);
z_edge=z_edge(find(z_edge(:,1)==Cycle(1) & abs(z_edge(:,2)-Cycle(2))<10^-6),:);
z_thick=z_thick(find(z_thick(:,1)==Cycle(1) & abs(z_thick(:,2)-Cycle(2))<10^-6),:);
load('Cycle Average Surface PAR.mat')
PAR_surf=SurfPAR(find(SurfPAR(:,1)==Cycle(1) & SurfPAR(:,2)==Cycle(2)),3);
clear surfPAR

wind = 5;   %m/s

%Setting the indices for state variables that sink or undergo mixing (mesozooplankton are assumed to be able to swim and hence maintain their position in the water column, they are not mixed)
fastsinkingIndices = [11,16,Thoriumon*(iThorium+9),N15on*(iN15+10)];
fastsinkingIndices = setdiff(fastsinkingIndices,0);
slowsinkingIndices = [10,15,Thoriumon*(iThorium+8),N15on*(iN15+9)];
slowsinkingIndices = setdiff(slowsinkingIndices,0);
if Thoriumon==1
    mesozooplanktonIndices = [4:7,26:29,35:38];
    lzdvmIndices = [5,27,36];
    pzdvmIndices = [7,29,38];
    livingIndices = [1:7,23:29,32:38];
else
    mesozooplanktonIndices = [4:7,25:28];
    lzdvmIndices = [5,27];
    pzdvmIndices = [7,29];
    livingIndices = [1:7,22:28];
end
mixedIndices = setdiff([1:length(tracer(1,:))],mesozooplanktonIndices);
nonlivingIndices = setdiff([1:length(tracer(1,:))],livingIndices);
lzdvmIndices = intersect(lzdvmIndices,[1:length(tracer(1,:))]);
pzdvmIndices = intersect(pzdvmIndices,[1:length(tracer(1,:))]);


%Setting boundary conditions
Deep=zeros(1,length(tracer(1,:)));
Deep(1,nonlivingIndices)=tracer(end,nonlivingIndices);
tracer(end,:)=[];

temperature(isnan(temperature(:,4)),:)=[];
if max(temperature(:,3))<200
    tmpdeep=min(temperature(:,4));
else
    tmpdeep=temperature(find(temperature(:,3)==200),4);
end
temperature=temperature(:,4);