clear all
close all


timelineplot=1;
movingvertplot=1;
movingoxygenplot=1;
movingvertratesplot=1;
CChlplot=1;
testconservation=0;
moving15Nplot=1;
movingDVMplot=1;

% Configuration
dt=1/24;				% time step (d)
dt_phys=1/24/10;         % time step for mixing and sinking (must be less than or equal to dt, and dt/dt_phys must be an integer)
dtz=1/24/100;           % time step for diel vertical migration (dt/dtz must be an integer)
sim_dur=30;				% simulation duration (days)


%PAR_surf=900;
Cycle = [0605,4];

CO2on=1;                  %Switch for turning on the CO2/Oxy submodule; set to 1 for on, 0 for off
Thoriumon=1;              %Switch for turning on the Thorium submodule; set to 1 for on, 0 for off
N15on=1;                  %Switch for turning on the N15 submodule; set to 1 for on, 0 for off
euponly=0;                %Switch for going from a euphotic zone only (DVM non-conserving) to a model that explicitly includes the mesopelagic and movement of vertical migrants

atmco2=400;
pH = 8;

%Initialization of state variables, etc.
[tracer,tracer_init,temperature,salinity,Deep,tmpdeep,Kz_mid,Kz_edge,z,z_edge,z_thick,PAR_surf,wind,iCO2,iThorium,iN15,fastsinkingIndices,slowsinkingIndices,mixedIndices,lzdvmIndices,pzdvmIndices,tot_iter]= ...
        InitializationDeep(Cycle,CO2on,Thoriumon,N15on,sim_dur,dt);
len = length(Kz_mid(:,1));
Kz_mid=[Kz_mid;[z(len+1:end,:),ones(size(z(len+1:end,1)))*10^-5*24*60*60]];  %Setting vertical eddy diffusivity beneath the 0.1% light level to 10^-5
len = length(Kz_edge(:,1));
Kz_edge=[Kz_edge;[z_edge(len+1:end,:),ones(size(z_edge(len+1:end,1)))*10^-5*24*60*60]];  %Setting vertical eddy diffusivity beneath the 0.1% light level to 10^-5

%Initializing parameters
[Param]=InitParameters(0);
omega_small=Param(95);
omega_large=Param(96);

%Calculating Coefficients for mixing and sinking
[MixingCoeff0,MixingCoeff1,MixingCoeff2,BottomCoeff]=CalculateMixingCoefficients(tracer(:,mixedIndices),z(:,3),z_edge(:,3),Kz_edge(:,4),dt_phys,Deep(:,mixedIndices));
[SinkingCoeff1,SinkingCoeff2,SinkingCoeff3]=CalculateSinkingCoefficients(tracer(:,[slowsinkingIndices,fastsinkingIndices]),z_edge(1:end-1,3),[omega_small*ones(size(slowsinkingIndices)),omega_large*ones(size(fastsinkingIndices))],dt_phys);



t=0
for i=1:tot_iter
    t = t+dt;
    if mod(t,1)>0.25 & mod(t,1)<0.75
        day=1;
    else
        day=0;
    end
    
    tic
    [tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality,basalexcretion_dvm,activexcretion_dvm,mortality_dvm]=NEMURObcp(tracer,Param,dt,PAR_surf,z(1:end-1,3),temperature(1:length(tracer(:,1))),salinity(1:length(tracer(:,1))),tmpdeep,day,CO2on,Thoriumon,N15on,euponly,iCO2,iThorium,iN15);
    tracer=tracer_out;
    if CO2on==1
        [tracer_out,co2flux,surfpco2,diffpco2] = Carbon(tracer,temperature(1),salinity(1),atmco2,z_thick(1,3),pH,wind,dt);
        tracer=tracer_out;
    end
    time1=toc;
    
    tic
    for j=1:dt/dt_phys
        [mixed_out]=mixing_ftcs(tracer(:,mixedIndices),Deep(:,mixedIndices),MixingCoeff0,MixingCoeff1,MixingCoeff2);   %No mixing for zooplankton
        tracer(:,mixedIndices)=mixed_out;
        [profile_new]=sinking_ftcs(tracer(:,[slowsinkingIndices,fastsinkingIndices]),SinkingCoeff1,SinkingCoeff2);
        tracer(:,[slowsinkingIndices,fastsinkingIndices])=profile_new;
    end
    time2=toc;
    
    
    if euponly==0
        if day==1
            targetdepth = z(min(find(PAR<10^-3)),3)+randn;   %10^-3 W / m2, from Bianchi et al. (2013)
        else
            phy=tracer(:,1)+tracer(:,2);
            targetdepth=sum(phy.*z(1:end-1,3))/sum(phy);   %Midpt of phytoplankton biomass
        end
        tic
        [DVMCoeff0,DVMCoeff1,DVMCoeff2]=CalculateDVMCoefficients(tracer(:,[lzdvmIndices,pzdvmIndices]),z(:,3),z_edge(:,3),dtz,targetdepth);
        time3=toc;
        tic
        for j=1:round(dt/dtz)
            [tracernew]=DVM_ftcs(tracer(:,[lzdvmIndices,pzdvmIndices]),DVMCoeff0,DVMCoeff1,DVMCoeff2);
            tracer(:,[lzdvmIndices,pzdvmIndices])=tracernew;
        end
        time4=toc;
    end
    
    %Graphing
    if timelineplot==1  PlotTimeline(t,tracer,z);   end
    if abs(round(t)-t)<dt/2
        if movingvertplot==1  PlotMovingVertical(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15);   end
        if testconservation==1    PlotConservation(tracer,tracer_init,z,t,N15on,iN15,Thoriumon,iThorium,euponly);   end
        if movingvertratesplot==1  PlotMovingVertRates(Cycle,tracer,GPP,NPP,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,mu_sp,mu_lp,m_chl,mu_chl,NO3up,NH4up,Siup,z,t,dt);   end
        if CChlplot==1    PlotCChl(Cycle,tracer,PAR,z,t,dt); end
        if moving15Nplot==1     PlotMoving15N(Cycle,tracer,z,t,iN15);   end
        if movingoxygenplot==1  PlotMovingOxygen(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15);   end
    end
    if movingDVMplot==1     PlotMovingDVM(Cycle,tracer,salinity,z,t,CO2on,Thoriumon,N15on,iCO2,iThorium,iN15,dt);   end
    
end




eupindex=max(find(PAR/PAR_surf>0.001));    %Index marking the 1% light level
eupdepth=z(eupindex,3);
ValidationRates=zeros(length(tracer(:,1)),14);
DVMRates=zeros(length(tracer(:,1)),2);
DayGraz=zeros(length(tracer(:,1)),2);
NightGraz=zeros(length(tracer(:,1)),2);
Sink=z_edge(2:end-1,3);
Sink(:,2:3)=0;
mixingloss=0;
sinkingloss=0;
CO2flux=0;
respiration_basal=0;
respiration_active=0;
mortality_DVM=0;
for i=1:1/dt
    t = t+dt;
    if mod(t,1)>0.25 & mod(t,1)<0.75
        day=1;
    else
        day=0;
    end
    [tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality,basalexcretion_dvm,activexcretion_dvm,mortality_dvm]=NEMURObcp(tracer,Param,dt,PAR_surf,z(1:end-1,3),temperature(1:length(tracer(:,1))),salinity(1:length(tracer(:,1))),tmpdeep,day,CO2on,Thoriumon,N15on,euponly,iCO2,iThorium,iN15);
    tracer=tracer_out;
    respiration_basal = respiration_basal + basalexcretion_dvm;
    respiration_active = respiration_active + activexcretion_dvm;
    mortality_DVM = mortality_DVM + mortality_dvm;
    
    if CO2on==1
        [tracer_out,co2flux,surfpco2,diffpco2] = Carbon(tracer,temperature(1),salinity(1),atmco2,z_thick(1,3),pH,wind,dt);
        tracer=tracer_out;
        CO2flux = CO2flux + co2flux;
    end
    
    for j=1:dt/dt_phys
        tracer_old = tracer;
        [mixed_out]=mixing_ftcs(tracer(:,mixedIndices),Deep(:,mixedIndices),MixingCoeff0,MixingCoeff1,MixingCoeff2);   %No mixing for zooplankton
        tracer(:,mixedIndices)=mixed_out;
        mixingloss = mixingloss + (sum(sum(tracer_old(1:eupindex,[1:3,10:13]))) - sum(sum(tracer(1:eupindex,[1:3,10:13]))));
        
        tracer_old = tracer;
        [profile_new]=sinking_ftcs(tracer(:,[slowsinkingIndices,fastsinkingIndices]),SinkingCoeff1,SinkingCoeff2);
        tracer(:,[slowsinkingIndices,fastsinkingIndices])=profile_new;
        sinkingloss = sinkingloss + (sum(sum(tracer_old(1:eupindex,[10:11]))) - sum(sum(tracer(1:eupindex,[10:11]))));
    end
    
    if euponly==0
        if day==1
            targetdepth = z(min(find(PAR<0.3*10^-3)),3);   %10^-3 W / m2, from Bianchi et al. (2013, GBC), modified because we use daily-average light which is lower (during the day) than instantaneous light used by Bianchi
            oxylim = z(min(find(tracer(:,19)<40)),3);        %40 mmol O m-3 seemed like a reasonable choice based on Supp. Fig. S7 of Bianchi et al. (2013, Nat. Geo.)
            targetdepth = min([targetdepth,oxylim])+randn;
        else
            phy=tracer(:,1)+tracer(:,2);
            targetdepth=sum(phy.*z(1:end-1,3))/sum(phy);   %Midpt of phytoplankton biomass
        end
        tic
        [DVMCoeff0,DVMCoeff1,DVMCoeff2]=CalculateDVMCoefficients(tracer(:,[lzdvmIndices,pzdvmIndices]),z(:,3),z_edge(:,3),dtz,targetdepth);
        time3=toc;
        tic
        for j=1:round(dt/dtz)
            [tracernew]=DVM_ftcs(tracer(:,[lzdvmIndices,pzdvmIndices]),DVMCoeff0,DVMCoeff1,DVMCoeff2);
            tracer(:,[lzdvmIndices,pzdvmIndices])=tracernew;
        end
        time4=toc;
    end
    
    if abs(mod(t,1)-0.5)<0.01
        PerPAR_noon=PAR/PAR_surf;
    end
    
    
    
end

tracer(end+1,:)=Deep;
