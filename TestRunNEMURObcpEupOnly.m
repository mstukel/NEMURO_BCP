clear all
close all


timelineplot=1;
movingvertplot=1;
movingvertratesplot=1;
CChlplot=1;
testconservation=9;

% Configuration
dt=1/24;				% time step (d)
dt_phys=1/24/10;         % time step for mixing and sinking (must be less than or equal to dt, and dt/dt_phys must be an integer
sim_dur=30;				% simulation duration (days)


%PAR_surf=900;
Cycle = [0605,5];

CO2on=1;                  %Switch for turning on the CO2/Oxy submodule; set to 1 for on, 0 for off
Thoriumon=1;              %Switch for turning on the Thorium submodule; set to 1 for on, 0 for off
N15on=1;                  %Switch for turning on the N15 submodule; set to 1 for on, 0 for off
euponly=1;                %Switch for going from a euphotic zone only (DVM non-conserving) to a model that explicitly includes the mesopelagic and movement of vertical migrants

atmco2=400;
pH = 8;

%Initialization of state variables, etc.
[tracer,tracer_init,temperature,salinity,Deep,tmpdeep,Kz_mid,Kz_edge,z,z_edge,z_thick,PAR_surf,wind,iCO2,iThorium,iN15,fastsinkingIndices,slowsinkingIndices,mixedIndices,tot_iter]= ...
        Initialization(Cycle,CO2on,Thoriumon,N15on,sim_dur,dt);


%Initializing parameters
[Param]=InitParameters(0);
omega_small=Param(95);
omega_large=Param(96);

% Param(42)=0
% Param(43)=0
% Param(44)=0
% Param(45:47)=0
% Param(51:53)=0
% Param(74:77)=0
% Param(82:84)=0
% Param(72)=0
% Param(73)=0
% Param(57)=0


[MixingCoeff0,MixingCoeff1,MixingCoeff2,BottomCoeff]=CalculateMixingCoefficients(tracer(:,mixedIndices),z(:,3),z_edge(:,3),Kz_edge(:,4),dt_phys,Deep(:,mixedIndices));
[SinkingCoeff1,SinkingCoeff2,SinkingCoeff3]=CalculateSinkingCoefficients(tracer(:,[slowsinkingIndices,fastsinkingIndices]),z_edge(1:end-1,3),[omega_small*ones(size(slowsinkingIndices)),omega_large*ones(size(fastsinkingIndices))],dt_phys);


tic
t=0
for i=1:tot_iter
    t = t+dt;
    if mod(t,1)>0.25 & mod(t,1)<0.75
        day=1;
    else
        day=0;
    end
    [tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality]=NEMURObcp(tracer,Param,dt,PAR_surf,z(1:end-1,3),temperature(1:length(tracer(:,1))),salinity(1:length(tracer(:,1))),tmpdeep,day,CO2on,Thoriumon,N15on,euponly,iCO2,iThorium,iN15);
    tracer=tracer_out;
    if CO2on==1
        [tracer_out,co2flux,surfpco2,diffpco2] = Carbon(tracer,temperature(1),salinity(1),atmco2,z_thick(1,3),pH,wind,dt);
        tracer=tracer_out;
    end
    
    for j=1:dt/dt_phys
        [mixed_out]=mixing_ftcs(tracer(:,mixedIndices),Deep(:,mixedIndices),MixingCoeff0,MixingCoeff1,MixingCoeff2);   %No mixing for zooplankton
        tracer(:,mixedIndices)=mixed_out;
        [profile_new]=sinking_ftcs(tracer(:,[slowsinkingIndices,fastsinkingIndices]),SinkingCoeff1,SinkingCoeff2);
        tracer(:,[slowsinkingIndices,fastsinkingIndices])=profile_new;
    end
    
    %Graphing
    if timelineplot==1  PlotTimeline(t,tracer,z);   end
    if abs(round(t)-t)<dt/2
        if movingvertplot==1  PlotMovingVertical(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15);   end
        if testconservation==1    PlotConservation(tracer,tracer_init,z,t,N15on,iN15,Thoriumon,iThorium,euponly);   end
        if movingvertratesplot==1  PlotMovingVertRates(Cycle,tracer,GPP,NPP,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,mu_sp,mu_lp,m_chl,mu_chl,NO3up,NH4up,Siup,z,t,dt);   end
        if CChlplot==1    PlotCChl(Cycle,tracer,PAR,z,t,dt); end
        
    end
    
end

time1=toc


tic
ValidationRates=zeros(length(tracer(:,1)),14);
DVMRates=zeros(length(tracer(:,1)),2);
DayGraz=zeros(length(tracer(:,1)),2);
NightGraz=zeros(length(tracer(:,1)),2);
Sink=z_edge(2:end-1,3);
Sink(:,2:3)=0;
mixingloss=0;
for i=1:1/dt
    t = t+dt;
    if mod(t,1)>0.25 & mod(t,1)<0.75
        day=1;
    else
        day=0;
    end
    [tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality]=NEMURObcp(tracer,Param,dt,PAR_surf,z(1:end-1,3),temperature(1:length(tracer(:,1))),salinity(1:length(tracer(:,1))),tmpdeep,day,CO2on,Thoriumon,N15on,euponly,iCO2,iThorium,iN15);
    tracer=tracer_out;
    ValidationRates=ValidationRates+[NPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,...
        LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,NO3up,NH4up,Siup];
    DVMRates=DVMRates+[DVM_excretion,DVM_mortality];
    if day==1
        DayGraz=DayGraz+[LZgrazchl,PZgrazchl];
    else
        NightGraz=NightGraz+[LZgrazchl,PZgrazchl];
    end
    if abs(mod(t,1)-0.5)<0.01
        PerPAR_noon=PAR/PAR_surf;
    end
    
    
    for i=1:dt/dt_phys
        
        [mixed_out,loss]=mixing(tracer(:,[1:3,8:16]),z(:,3),z_edge(:,3),Kz_edge(:,4),dt_phys,Deep(:,[1:3,8:16]));   %No mixing for zooplankton
        mixingloss=mixingloss+sum(-loss(:,[1:3,6:9]));   %Only organic matter containing groups - this is at the bottom of the model
        tracer(:,[1:3,8:16])=mixed_out;
        [profile_new,Flux]=sinking(tracer(:,10),z_edge(1:end-1,3),omega_small,dt_phys);
        Sink(:,2)=Sink(:,2)+Flux;
        tracer(:,10)=profile_new;
        [profile_new,Flux]=sinking(tracer(:,11),z_edge(1:end-1,3),omega_large,dt_phys);
        Sink(:,2)=Sink(:,2)+Flux;
        tracer(:,11)=profile_new;
        [profile_new,Flux]=sinking(tracer(:,15),z_edge(1:end-1,3),omega_small,dt_phys);
        Sink(:,3)=Sink(:,3)+Flux;
        tracer(:,15)=profile_new;
        [profile_new,Flux]=sinking(tracer(:,16),z_edge(1:end-1,3),omega_large,dt_phys);
        Sink(:,3)=Sink(:,3)+Flux;
        tracer(:,16)=profile_new;
    end
end
time2=toc
tracer(end+1,:)=Deep;

mixingloss