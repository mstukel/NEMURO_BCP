function [tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality,basalexcretion_dvm,activexcretion_dvm,mortality_dvm]=NEMURObcp(tracer,Param,dt,PAR_surf,depth_t,temperature,salinity,tmpdeep,day,CO2on,Thoriumon,N15on,euponly,iCO2,iThorium,iN15)
%[tracer_out,NPP,GPP,mu_chl,m_chl,mu_sp,m_sp,mu_lp,m_lp,SP2SZ,LP2SZ,LZresgraz,LZdvmgraz,PZresgraz,PZdvmgraz,LZgrazchl,PZgrazchl,NO3up,NH4up,Siup,PAR,DVM_excretion,DVM_mortality]=NEMURObcp(tracer,Param,dt,PAR_surf,depth_t,temperature,salinity,tmpdeep,day,CO2on,Thoriumon,N15on,euponly)

% *************************************************************************** %
% North Pacific Ecosystem Model for Understanding Regional Oceanography - Biological Carbon Pump (NEMURObcp)
% Description: 1-D Lower Trophic Level Marine Biogeochemical Numerical Model modified to investigate Biological Pump
% Modified from NEMURO (Kishi et al. 2007) and NEMURO-GoM (Shropshire et al. 2020)
% This is written in a backward Euler implementation to ensure stability
% Written By: Mike Stukel
% *************************************************************************** %


% ********** Parameters ********** %

% Light Attenuation
par_frac=1.0;                  % fraction of incident radiation that is photosynthetically active radiation
ext_w=0.03;                     % attenuation due to water - [0.04]
ext_chl = Param(101);                  % attenuation due to phytoplankton chl - [0.05]


% Temperature Dependence
tdep1=Param(1); scale_tmp_phyto1=1;        % phytoplankton growth
tdep2=Param(1); scale_tmp_phyto2=1;        % phytoplankton mortality
tdep3=Param(1); scale_tmp_phyto3=1;        % phytolpankton respiration
tdep4=Param(1); scale_tmp_zoo1=1;          % zooplankton grazing (1:0.0693, 0.14:0.1386, 0.016:0.2079)
tdep5=Param(1); scale_tmp_zoo2=1;          % zooplankton mortality (1:0.0693, 0.14:0.1386)
tdep6=Param(1); scale_tmp_decomp=1;        % decomposition (1:0.0693, 0.14:0.1386)

% Small Phytoplankton
vmax_sp = Param(2);                  % max growth rate at 0 deg C (1/d) - [0.4]
k_NO_sp = Param(3);                  % nitrate half saturation constant - [1.0]
k_NH_sp = Param(4);                  % ammonium half saturation constant - [0.1]
alpha_sp = Param(5);    	    	% PI curve parameter alpha - [0.01]
beta_sp = Param(6);                 % PI curve parmater beta - [4.5e-4], 1.4e-3
ref_resp_sp = Param(7);            % respiration at 0 deg C (1/d) - [0.03]
ref_mort_sp = Param(8);              % mortality at at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
ext_excr_sp = Param(9);             % extracellular excretion - [0.135]
tc_v_sp = Param(1);                % tmp coefficent for growth - [0.0693]
tc_r_sp = Param(1);                % tmp coefficent for respiration - [0.0519]
tc_m_sp = Param(1);                % tmp coeffiect for mortality - [0.0693]
inh_NH_NO_sp = Param(10);               % nitrate uptake inhibiton, NOTE: parameter currently not being used - [1.4]

% Large Phytoplankton
vmax_lp = Param(11);                  % max growth rate at 0 deg C (1/d) - [0.8]
k_NO_lp = Param(12);                  % nitrate half saturation constant - [3.0]
k_NH_lp = Param(13);                  % ammonium half saturation constant - [0.3]
k_SI_lp = Param(14);                  % silica half saturation constant - [6.0]
alpha_lp = Param(15);                 % PI curve parameter alpha - [0.01], 1.4e-3
beta_lp= Param(16);                  % PI curve parmater beta - [4.5e-4]
ref_resp_lp = Param(17);             % respiration at 0 deg C (1/d) - [0.03]
ref_mort_lp = Param(18);              % mortality at at 0 deg C (1/d) - [0.029], Linear Mort - ~[0.001]
ext_excr_lp = Param(19);              % extracellular excretion - [0.0135]
tc_v_lp = Param(1);                % tmp coefficent for growth - [0.0693]
tc_r_lp = Param(1);                % tmp coefficent for respiration - [0.0519]
tc_m_lp = Param(1);                % tmp coeffiect for mortality - [0.0693]
inh_NH_NO_lp = Param(20);               % nitrate uptake inhibiton, NOTE: parameter currently not being used - [1.4]


% Small Zooplankton
gmax_sz_sp = Param(21);               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.04]
gmax_sz_lp = Param(22);              % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.04]
tc_g_sz_sp = Param(1);             % tmp coefficent for grazing - [0.0693]
tc_g_sz_lp = Param(1);             % tmp coefficent for grazing - [0.0693]
iv_sz_sp = Param(23);                 % ivlev constant - [1.4]
iv_sz_lp = Param(24);                 % ivlev constant - [1.4]
thresh_sz_sp = Param(25);           % feeding threshold on small phytoplankton - [0.043]
thresh_sz_lp = Param(26);           % feeding threshold on small phytoplankton - [0.043]
ref_mort_sz = Param(27);              % mortality at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
tc_m_sz = Param(1);                % tmp coefficent for mortality - [0.0693]
ae_sz = Param(28);                   % assimilation efficency - [0.7]
gge_sz = Param(29);                  % gross growth efficency - [0.3]

% Large Zooplankton - Resident
gmax_lzres_sp = Param(30);           % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.1]
gmax_lzres_lp = Param(31);               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.4]
gmax_lzres_sz = Param(32);               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.4]
tc_g_lzres_sp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_lzres_lp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_lzres_sz = Param(1);             % tmp coefficent for grazing -[0.0693]
iv_lzres_sp = Param(33);                      % ivlev constant - [1.4]
iv_lzres_lp = Param(34);                 % ivlev constant - [1.4]
iv_lzres_sz = Param(35);                 % ivlev constant - [1.4]
ae_lzres = Param(36);                   % assimilation efficency - [0.7]
act_res_lzres = Param(37);                  % active respiration (fraction of grazing) - [0.3]
mort_day_lzres = Param(38);                           % Quadratic Mortality term for LZ during the day
thresh_lzres_sp = Param(39);                  % feeding threshold on small phytoplankton - [0.04]
thresh_lzres_lp = Param(40);            % feeding threshold on large phytoplankton - [0.04]
thresh_lzres_sz = Param(41);            % feeding threshold on small zooplankton  - [0.04]
mort_night_lzres = Param(42);           % Quadratic Mortality term for LZ during the night
tc_m_lz = Param(1);                % tmp coefficent for mortality - [0.0693]
Ikeda_a2 = Param(43);
Ikeda_lz = Param(44);

% Large Zooplankton   -   DVM
gmax_lzdvm_sp = Param(45);           % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.1]
gmax_lzdvm_lp = Param(46);               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.4]
gmax_lzdvm_sz = Param(47);               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.4]
tc_g_lzdvm_sp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_lzdvm_lp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_lzdvm_sz = Param(1);             % tmp coefficent for grazing -[0.0693]
iv_lzdvm_sp = Param(48);                        % ivlev constant - [1.4]
iv_lzdvm_lp = Param(49);                 % ivlev constant - [1.4]
iv_lzdvm_sz = Param(50);                 % ivlev constant - [1.4]
ae_lzdvm = Param(51);                   % assimilation efficency - [0.7]
act_res_lzdvm = Param(52);                  % active respiration (fraction of grazing) - [0.3]
mort_day_lzdvm = Param(53);                           % Quadratic Mortality term for LZ during the night
thresh_lzdvm_sp = Param(39);                  % feeding threshold on small phytoplankton - [0.04]
thresh_lzdvm_lp = Param(40);            % feeding threshold on large phytoplankton - [0.04]
thresh_lzdvm_sz = Param(41);            % feeding threshold on small zooplankton  - [0.04]
mort_night_lzdvm = Param(42);          % Quadratic Mortality term for LZ during the night

% Predatory Zooplankton   -   Resident
gmax_pzres_sp = Param(54);               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.2]
gmax_pzres_lp = Param(55);               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.2]
gmax_pzres_sz = Param(56);               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.2]
gmax_pzres_lz = Param(57);               % max grazing rate on large zooplankton at 0 deg C (1/d) - [0.2]
tc_g_pzres_sp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzres_lp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzres_sz = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzres_lz = Param(1);             % tmp coefficent for grazing -[0.0693]
iv_pzres_sp = Param(58);                 % ivlev constant - [1.4]
iv_pzres_lp = Param(59);                 % ivlev constant - [1.4]
iv_pzres_sz = Param(60);                 % ivlev constant - [1.4]
iv_pzres_lz = Param(61);                 % ivlev constant - [1.4]
ae_pzres = Param(62);                    % assimilation efficency - [0.7]
act_res_pzres = Param(63);                  % active respiration (fraction of grazing)  - [0.3]
mort_day_pzres = Param(64);            % Quadratic Mortality term for LZ during the night
inh_szlzlp_sp = Param(65);            % grazing on large phytoplankton inhibition by small and large zooplankton - [4.605]
inh_szlz_lp = Param(66);            % grazing on large phytoplankton inhibition by small and large zooplankton - [4.605]
inh_lz_sz = Param(67);               % grazing on small zooplankton inhibition by large zooplankton - [3.01]
thresh_pzres_sp = Param(68);            % feeding threshold on small phytoplankton - [0.04]
thresh_pzres_lp = Param(69);            % feeding threshold on large phytoplankton - [0.04]
thresh_pzres_sz = Param(70);            % feeding threshold on small zooplankton - [0.04]
thresh_pzres_lz = Param(71);            % feeding threshold on large zooplankton - [0.04]
mort_night_pzres = Param(72);             % Quadratic Mortality term for LZ during the day
Ikeda_pz = Param(73);
tc_m_pz = Param(1);                % tmp coefficent for mortality - [0.0693]



% Predatory Zooplankton   -   DVM
gmax_pzdvm_sp = Param(74);               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.2]
gmax_pzdvm_lp = Param(75);               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.2]
gmax_pzdvm_sz = Param(76);               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.2]
gmax_pzdvm_lz = Param(77);               % max grazing rate on large zooplankton at 0 deg C (1/d) - [0.2]
tc_g_pzdvm_sp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzdvm_lp = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzdvm_sz = Param(1);             % tmp coefficent for grazing -[0.0693]
tc_g_pzdvm_lz = Param(1);             % tmp coefficent for grazing -[0.0693]
iv_pzdvm_sp = Param(78);                 % ivlev constant - [1.4]
iv_pzdvm_lp = Param(79);                 % ivlev constant - [1.4]
iv_pzdvm_sz = Param(80);                 % ivlev constant - [1.4]
iv_pzdvm_lz = Param(81);                 % ivlev constant - [1.4]
ae_pzdvm = Param(82);                    % assimilation efficency - [0.7]
act_res_pzdvm = Param(83);                  % active respiration (fraction of grazing)  - [0.3]
mort_day_pzdvm = Param(84);                           % Quadratic Mortality term for LZ during the day
thresh_pzdvm_sp = Param(68);                      % feeding threshold on small phytoplankton - [0.04]
thresh_pzdvm_lp = Param(69);            % feeding threshold on large phytoplankton - [0.04]
thresh_pzdvm_sz = Param(70);            % feeding threshold on small zooplankton - [0.04]
thresh_pzdvm_lz = Param(71);            % feeding threshold on large zooplankton - [0.04]
mort_night_pzdvm = Param(72);                           % Quadratic Mortality term for LZ during the night



% Nutrients (11)
ref_nitr = Param(85);                % nitrification at 0 deg C (1/d) - [0.03]
ref_dec_PON_NH = Param(86);           % decompositon from PON to NH at 0 deg C (1/d) - [0.1]
ref_dec_PON_DON = Param(87);          % decompositon from PON to DON at 0 deg C (1/d) - [0.1]
ref_dec_PON_LPON = Param(88);          % aggregation PON to LON at 0 deg C (1/d) - [0.1]
ref_dec_LPON_NH = Param(89);           % decompositon from LPON to NH at 0 deg C (1/d) - [0.1]
ref_dec_LPON_DON = Param(90);          % decompositon from LPON to DON at 0 deg C (1/d) - [0.1]
ref_dec_DON_NH = Param(91);           % decompositon from DON to NH at 0 deg C (1/d) - [0.2]
f_DON_DONref = Param(92);           % fraction of labile DON that gets converted to refractory DON
ref_dec_DONref_NH = Param(93);           % decompositon from DON to NH at 0 deg C (1/d) - [0.2]
ref_dec_OP_SI = Param(94);            % decompositon from OP to SI at 0 deg C (1/d) - [0.1]
tc_nitr = Param(1);                % tmp coefficent for nitrification - [0.0693]
tc_dec_PON_NH = Param(1);          % tmp coefficent for decomposition  - [0.0693]
tc_dec_PON_DON = Param(1);         % tmp coefficent for decomposition - [0.0693]
tc_dec_PON_LPON = Param(1);         % tmp coefficent for decomposition - [0.0693]
tc_dec_DON_NH = Param(1);          % tmp coefficent for decomposition - [0.0693]
tc_dec_OP_SI = Param(1);           % tmp coefficent for decomposition - [0.0693]
r_SI_N = 1.0;                   % ratio of silica to nitrogen (SI/N) - [2.0]
sink=Param(95);                    		% sinking speed per day - [40]
Lsink=Param(96);                    		% sinking speed per day for large PON and OP - [40]

koxy_mic = Param(97);               %Oxygen limitation half saturation constant (microbes = heterotrophic protists + implicit bacteria)
koxy_met = Param(98);               %Oxygen limitation half saturation constant (metazoans)

% Chl2c Submodel (Li et al., 2010)
chl2c_sp_min = Param(99); 		% minimum Chl:C ratio - [0.000]
chl2c_lp_min = Param(100); 		% maximum Chl:C ratio - [0.005] 
chl2c_sp_max = Param(101); 		% minimum Chl:C ratio - [0.03]
chl2c_lp_max = Param(102); 		% maximum Chl:C ratio - [0.061] 

%Constant ratios
carbon_nitrogen = 106/16;                 % Redfield C:N ratio (mol:mol)
oxy_nitrogen_nh4 = 118/16;     % Redfield O:N ratio for NH4 assimilation or regeneration
oxy_nitrogen_no3 = 150/16;     % Redfield O:N ratio for NO3 assimilation
nitrogen_phosphorus = 16;  %Redfield N:P ratio


% ######################################################## %
% #################### NEMURObcp ######################### %
% ######################################################## %





%Set variables from tracer array
%_v1 is prior value of variable
%_v2 is updated value of variable
sp_v1=tracer(:,1); sp_v2=sp_v1;                 %SP = non-diatom (small) phytoplankton
lp_v1=tracer(:,2); lp_v2=lp_v1;                 %LP = Diatoms
sz_v1=tracer(:,3); sz_v2=sz_v1;                 %SZ = small zooplankton (protistan zooplankton)
lzres_v1=tracer(:,4); lzres_v2=lzres_v1;        %LZRES = resident large zooplankton (<1-mm mesozooplankton that do not vertically migrate)
lzdvm_v1=tracer(:,5); lzdvm_v2=lzdvm_v1;        %LZDVM = Vertically-migrating large zooplankton (<1-mm mesozooplankton that vertically migrate every day)
pzres_v1=tracer(:,6); pzres_v2=pzres_v1;        %PZRES = resident "predatory" zooplankton (>1-mm mesozooplankton that do not vertically migrate)
pzdvm_v1=tracer(:,7); pzdvm_v2=pzdvm_v1;        %PZDVM = Vertically-migrating "predatory" zooplankton (?1-mm mesozooplankton that vertically migrate every day)
NO_v1=tracer(:,8); NO_v2=NO_v1;                 %NO = Nitrate
NH_v1=tracer(:,9); NH_v2=NH_v1;                 %NH = Ammonium
PON_v1=tracer(:,10); PON_v2=PON_v1;             %PON = Particulate organic nitrogen (small non-living detritus)
LPON_v1=tracer(:,11); LPON_v2=LPON_v1;          %LPON = Large particulate organic nitrogen (non-living detritus with a higher sinking speed; e.g., fecal pellets and aggregates)
DON_v1=tracer(:,12); DON_v2=DON_v1;             %DON = Dissolved organic nitrogen (labile)
DONref_v1=tracer(:,13); DONref_v2=DONref_v1;    %DONref = Dissolved organic nitrogen (refractory)
SI_v1=tracer(:,14); SI_v2=SI_v1;                %SI = silicic acid
OP_v1=tracer(:,15); OP_v2=OP_v1;                %OP = Opal (small non-living biogenic silica; i.e., non-living diatom frustules)
LOP_v1=tracer(:,16); LOP_v2=LOP_v1;             %LOP = Large opal (large non-living biogenic silica; e.g., silica in fecal pellets and aggregates)
Chl_ps_v1=tracer(:,17); Chl_ps_v2=Chl_ps_v1;    %Chlorophyll associated with small phytoplankton (units are mg Chl a m^-3)
Chl_pl_v1=tracer(:,18); Chl_pl_v2=Chl_pl_v1;    %Chlorophyll associated with large phytoplankton (units are mg Chl a m^-3)
oxy_v1=tracer(:,19); oxy_v2=oxy_v1;             %Oxygen
if CO2on==1
    DIC_v1=tracer(:,iCO2+0); DIC_v2=DIC_v1;             %Dissolved inorganic carbon
    alk_v1=tracer(:,iCO2+1); alk_v2=alk_v1;             %Alkalinity
end

    
tmp=temperature(:,1);



%%%%%%%%%% Backward Euler Loop %%%%%%%%%%%%%%
for iter=1:3
    
    % ---- Determine Temperature Light ---- %
    
    for zz=1:length(tracer(:,1))       %Looping over all depth layers
        if zz==1
            rad=PAR_surf*par_frac*exp(-1*(ext_chl*(Chl_ps_v2(zz)+Chl_pl_v2(zz))+ext_w)*depth_t(zz));
        else
            rad=PAR(zz-1,1)*exp(-1*(ext_chl*(Chl_ps_v2(zz)+Chl_pl_v2(zz))+ext_w)*(depth_t(zz)-depth_t(zz-1)));
        end
        PAR(zz,1)=rad;
    end
    rad=PAR;
    
    %****** -------- Chl:C Sub Model ----------------------- ******* %
    
    alpha_chl_sp = alpha_sp./chl2c_sp_max; 		% chlorophyll specific initial slope of Pâ€?E curve - [0.28 +- 0.11]
    alpha_chl_lp = alpha_lp./chl2c_lp_max; 		% chlorophyll specific initial slope of Pâ€?E curve - [0.28 +- 0.11]
    
    %Small Phytoplankton
    t_lim_v_sp = exp(tc_v_sp*tmp);
    NO_lim_sp = NO_v2./(NO_v2+k_NO_sp).*exp(-inh_NH_NO_sp*NH_v2);
    NH_lim_sp = NH_v2./(NH_v2+k_NH_sp);
    PmC = vmax_sp*t_lim_v_sp.*(NO_lim_sp+NH_lim_sp).*(alpha_sp/(alpha_sp+beta_sp)).*(beta_sp/(alpha_sp+beta_sp)).^(beta_sp/alpha_sp);
    omegageider = chl2c_sp_max.*alpha_chl_sp./PmC.*rad;
    chl2c_sp= chl2c_sp_max./(1 + 0.5*omegageider);
    chl2c_sp=max(chl2c_sp,chl2c_sp_min);
    Chl_ps_v2 = sp_v2 .* carbon_nitrogen .* 12 .* chl2c_sp;
    
    %Large Phytoplankton
    t_lim_v_lp = exp(tc_v_lp*tmp);
    NO_lim_lp = NO_v2./(NO_v2+k_NO_lp).*exp(-inh_NH_NO_lp*NH_v2);
    NH_lim_lp = NH_v2./(NH_v2+k_NH_lp);
    SI_lim_lp = SI_v2./(SI_v2+k_SI_lp);
    NS_lim_comp = min(1.0,SI_lim_lp./(NO_lim_lp+NH_lim_lp));
    PmC=vmax_lp*t_lim_v_lp.*(NO_lim_lp+NH_lim_lp).*NS_lim_comp.*(alpha_lp/(alpha_lp+beta_lp))*(beta_lp/(alpha_lp+beta_lp))^(beta_lp/alpha_lp);
    omegageider = chl2c_lp_max.*alpha_chl_lp./PmC.*rad;
    chl2c_lp= chl2c_lp_max./(1 + 0.5*omegageider);
    chl2c_lp=max(chl2c_lp,chl2c_lp_min);
    Chl_pl_v2 = lp_v2 * carbon_nitrogen * 12 .* chl2c_lp;
    
    
    % ### Functional Groups: ### %
    
    % ************ Nutrient Uptake Terms ****************** %
    
    %-----------------------------%
    % --- Small Phytoplankton --- %
    %-----------------------------%
    
    % Calculate Limitations
    t_lim_v_sp = exp(tc_v_sp*tmp)*scale_tmp_phyto1;
    l_lim_sp = (1.0-exp(-1.0*alpha_sp*rad/vmax_sp)).*exp(-1.0*beta_sp.*rad/vmax_sp);
    NO_lim_sp = NO_v2./(NO_v2+k_NO_sp).*exp(-inh_NH_NO_sp*NH_v2);
    NH_lim_sp = NH_v2./(NH_v2+k_NH_sp);
    
    % sp NO Uptake
    NO_v2 = NO_v1./(1.0+dt./NO_v2.*(vmax_sp*t_lim_v_sp...
        .*l_lim_sp.*NO_lim_sp.*sp_v2));
    delta_NO = NO_v1-NO_v2;
    
    % sp NH Uptake
    NH_v2 = NH_v1./(1.0+dt./NH_v2.*(vmax_sp*t_lim_v_sp...
        .*l_lim_sp.*NH_lim_sp.*sp_v2));
    delta_NH = NH_v1 - NH_v2;
    
    % sp Excretion
    delta_sp = ext_excr_sp*(delta_NO+delta_NH);
    
    % Update
    DON_v2 = DON_v1+delta_sp;
    sp_v2 = sp_v1+delta_NO+delta_NH-delta_sp;
    
    %Calculate f-ratio for SP
    np_frac_sp = delta_NO./(delta_NO+delta_NH);
    np_frac_sp(isnan(np_frac_sp))=0;  %Note this step is necessary for low light regions where delta_NO and delta_NH can both be equal to 0
    
    % Diagnostic
    gpp_sp=delta_NO+delta_NH;
    npp_sp=delta_NO+delta_NH-delta_sp;
    mu_sp=npp_sp./sp_v1;                %Note that mu_sp(zz,1) gets modified later in the code
    NO3up=delta_NO;                     %Note that NO3up(zz,1) gets modified later in the code
    NH4up=delta_NH;                     %Note that NH4up(zz,1) gets modified later in the code
    NO3up_sp=delta_NO;
    NH4up_sp=delta_NH;
    exu_sp=delta_sp;
    % ---------------------- %
    
    %----------------------------%
    % --- Large Phytoplankton ---%
    %----------------------------%
    
    % Calculate Limitations
    t_lim_v_lp = exp(tc_v_lp*tmp)*scale_tmp_phyto1;
    l_lim_lp = (1.0-exp(-1.0*alpha_lp*rad/vmax_lp)).*exp(-1.0*beta_lp*rad/vmax_lp);
    NO_lim_lp = NO_v2./(NO_v2+k_NO_lp).*exp(-inh_NH_NO_lp*NH_v2);
    NH_lim_lp = NH_v2./(NH_v2+k_NH_lp);
    SI_lim_lp = SI_v2./(SI_v2+k_SI_lp);
    NS_lim_comp = min(1.0,SI_lim_lp./(NO_lim_lp+NH_lim_lp));
    
    % lp NO Uptake
    cff_NO = NO_v2-(NO_v2./(1.0+dt./NO_v2.*(vmax_lp*t_lim_v_lp.*l_lim_lp...
        .*NO_lim_lp.*NS_lim_comp.*lp_v2)));
    
    % lp NH Uptake
    cff_NH = NH_v2-(NH_v2./(1.0+dt./NH_v2.*(vmax_lp*t_lim_v_lp.*l_lim_lp...
        .*NH_lim_lp.*NS_lim_comp.*lp_v2)));
    
    cff_SI = SI_v2-(SI_v2./(1.0+dt./SI_v2.*...
        (r_SI_N.*(vmax_lp.*t_lim_v_lp.*l_lim_lp.*NH_lim_lp.*NS_lim_comp.*lp_v2)+...
        r_SI_N.*(vmax_lp.*t_lim_v_lp.*l_lim_lp.*NO_lim_lp.*NS_lim_comp.*lp_v2))     ));
    
    delta_NO = min(cff_NO,cff_NO.*cff_SI./(cff_NH+cff_NO));
    delta_NH = min(cff_NH,cff_NH.*cff_SI./(cff_NH+cff_NO));
    
    % lp Excretion
    delta_lp = ext_excr_lp*(delta_NO+delta_NH);
    
    % Update
    NO_v2 = NO_v2-delta_NO;
    NH_v2 = NH_v2-delta_NH;
    SI_v2 = SI_v1-r_SI_N.*(delta_NO+delta_NH)+r_SI_N*delta_lp;
    DON_v2 = DON_v2+delta_lp;
    lp_v2 = lp_v1+delta_NO+delta_NH-delta_lp;
    
    %Calculate f-ratio for LP
    np_frac_lp = delta_NO./(delta_NO+delta_NH);
    np_frac_lp(isnan(np_frac_lp))=0;  %Note this step is necessary for low light regions where delta_NO and delta_NH can both be equal to 0
    
    % Diagnostic
    gpp_lp=delta_NO+delta_NH;
    npp_lp=delta_NO+delta_NH-delta_lp;
    mu_lp=npp_lp./lp_v1;                %Note that mu_sp(zz,1) gets modified later in the code
    NO3up=NO3up+delta_NO;
    NH4up=NH4up+delta_NH;
    Siup=r_SI_N*(delta_NO+delta_NH)-r_SI_N*delta_lp;
    NO3up_lp=delta_NO;
    NH4up_lp=delta_NH;
    exu_lp=delta_lp;
    
    % ---------------------- %
    % *********** END Nutrient Uptake Terms ********** %
    
    %------------------------------------------------------%
    % ---- Phytoplankton Respiration & Mortality Terms --- %
    %------------------------------------------------------%
    
    % sp respiration
    t_lim_r_sp = exp(tc_r_sp*tmp)*scale_tmp_phyto3;
    respiration_sp = sp_v2-(sp_v2./(1.0+dt./sp_v2.*(ref_resp_sp*t_lim_r_sp.*sp_v2)));
    
    % lp respiration
    t_lim_r_lp = exp(tc_r_lp*tmp)*scale_tmp_phyto3;
    respiration_lp = lp_v2-(lp_v2./(1.0+dt./lp_v2.*(ref_resp_lp*t_lim_r_lp.*lp_v2)));
    
    % Update
    np_frac_sp = np_frac_sp.*min(gpp_sp./respiration_sp,1);   %This code fixes the glitch in NEMURO of phytoplankton producing nitrate when respiration > primary production
    sp_v2 = sp_v2-respiration_sp;
    NO_v2 = NO_v2+respiration_sp.*np_frac_sp;
    NH_v2 = NH_v2+respiration_sp.*(1.0-np_frac_sp);
    mu_sp=mu_sp-respiration_sp./sp_v1;
    
    np_frac_lp = np_frac_lp.*min(gpp_sp./respiration_sp,1);   %This code fixes the glitch iNEMURO of phytoplankton producing nitrate when respiration > primary production
    lp_v2 = lp_v2-respiration_lp;
    NO_v2 = NO_v2+respiration_lp.*np_frac_lp;
    NH_v2 = NH_v2+respiration_lp.*(1.0-np_frac_lp);
    SI_v2 = SI_v2+respiration_lp*r_SI_N;
    mu_lp=mu_lp-respiration_lp./lp_v1;
    
    %Diagnostic
    npp_lp = npp_lp - respiration_lp;
    npp_sp = npp_sp - respiration_sp;
    NPP = npp_lp + npp_sp;
    GPP = gpp_lp + gpp_sp;
    
    % sp Mortality
    t_lim_m_sp = exp(tc_m_sp*tmp)*scale_tmp_phyto2;
    mortality_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(ref_mort_sp*t_lim_m_sp.*sp_v2)));
    
    % lp Mortality
    t_lim_m_lp = exp(tc_m_lp*tmp)*scale_tmp_phyto2;
    mortality_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(ref_mort_lp*t_lim_m_lp.*lp_v2)));
    
    % Update
    sp_v2 = sp_v2-mortality_sp;
    lp_v2 = lp_v2-mortality_lp;
    PON_v2 = PON_v1+mortality_sp+mortality_lp;
    OP_v2 = OP_v1+mortality_lp*r_SI_N;
    %mu_sp=mu_sp-mortality_sp./sp_v1;
    %mu_lp=mu_lp-mortality_lp./lp_v1;
    mu_chl=(mu_sp.*Chl_ps_v2 + mu_lp.*Chl_pl_v2)./(Chl_ps_v2+Chl_pl_v2);
    
    % ------------------------- %
    
    
    % ********* Zooplankton Grazing Terms **********  %
    
    %-------------------------------------------%
    % ------- Small Zooplankton Grazing ------- %
    %-------------------------------------------%
    
    % sz grazing on sp
    t_lim_g_sz_sp = exp(tc_g_sz_sp*tmp)*scale_tmp_zoo1;
    g_lim_sz_sp = max(0,1-exp(iv_sz_sp*(thresh_sz_sp-sp_v2)));
    oxy_lim_mic = oxy_v2./(koxy_mic+oxy_v2);
    delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_sz_sp*t_lim_g_sz_sp.*g_lim_sz_sp.*oxy_lim_mic.*sz_v2)));
    
    % sz grazing on lp
    t_lim_g_sz_lp = exp(tc_g_sz_lp*tmp)*scale_tmp_zoo1;
    g_lim_sz_lp = max(0,1-exp(iv_sz_lp*(thresh_sz_lp-lp_v2)));
    delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_sz_lp*t_lim_g_sz_lp.*g_lim_sz_lp.*oxy_lim_mic.*sz_v2)));
    
    % Update
    sp_v2 = sp_v2-delta_sp;
    lp_v2 = lp_v2-delta_lp;
    sz_v2 = sz_v1+gge_sz*(delta_sp+delta_lp);
    NH_v2 = NH_v2+(ae_sz-gge_sz)*(delta_sp+delta_lp);
    PON_v2 = PON_v2+(1-ae_sz)*(delta_sp+delta_lp);
    OP_v2 = OP_v2+delta_lp*r_SI_N;
    
    % Diagnostics
    SP2SZ=delta_sp;
    LP2SZ=delta_lp;
    m_sp = delta_sp./sp_v1;   %Specific mortality rate of SP due to protistan (SZ) grazing
    m_lp = delta_lp./lp_v1;   %Specific mortality rate of LP due to protistan (SZ) grazing
    m_chl=(m_sp.*Chl_ps_v2 + m_lp.*Chl_pl_v2)./(Chl_ps_v2+Chl_pl_v2);
    % ------------------------%
    
    
    if day==1 & euponly==1      %Mesozooplankton grazing terms during the day
        
        %----------------------------------------------%
        % --- Large Zooplankton (Resident) Grazing --- %
        %----------------------------------------------%
        
        % lzres grazing on sp
        t_lim_g_lzres_sp =exp(tc_g_lzres_sp*tmp)*scale_tmp_zoo1;
        g_lim_lzres_sp = max(0,1-exp(iv_lzres_sp*(thresh_lzres_sp-sp_v2)));
        oxy_lim_met = oxy_v2./(koxy_met+oxy_v2);
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_lzres_sp*t_lim_g_lzres_sp.*g_lim_lzres_sp.*oxy_lim_met.*lzres_v2)));
        
        % lzres grazing on lp
        t_lim_g_lzres_lp =exp(tc_g_lzres_lp*tmp)*scale_tmp_zoo1;
        g_lim_lzres_lp = max(0,1-exp(iv_lzres_lp*(thresh_lzres_lp-lp_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_lzres_lp*t_lim_g_lzres_lp.*g_lim_lzres_lp.*oxy_lim_met.*lzres_v2)));
        
        % lzres grazing on sz
        t_lim_g_lzres_sz = exp(tc_g_lzres_sz*tmp)*scale_tmp_zoo1;
        g_lim_lzres_sz = max(0,1-exp(iv_lzres_sz*(thresh_lzres_sz-sz_v2)));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_lzres_sz*t_lim_g_lzres_sz.*g_lim_lzres_sz.*oxy_lim_met.*lzres_v2)));
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzres_v2 = lzres_v1+(ae_lzres-act_res_lzres)*(delta_sp+delta_lp+delta_sz);
        NH_v2 = NH_v2+(act_res_lzres)*(delta_sp+delta_lp+delta_sz);
        LPON_v2 = LPON_v1+(1-ae_lzres)*(delta_sp+delta_lp+delta_sz);
        LOP_v2 = LOP_v1+delta_lp*r_SI_N;
        
        
        % Diagnostics
        SP2LZres=delta_sp;
        LP2LZres=delta_lp;
        SZ2LZres=delta_sz;
        SP2LZdvm=0;
        LP2LZdvm=0;
        SZ2LZdvm=0;
        LZresgraz=delta_sp+delta_lp;
        LZdvmgraz=zeros(size(lzdvm_v2));
        LZgrazchl = delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;
        % -------------------------%
        
        
        %----------------------------------------------%
        % ------- Predatory Zooplankton Grazing ------ %
        %----------------------------------------------%
        
        % pzres grazing on sp
        t_lim_g_pzres_sp = exp(tc_g_pzres_sp*tmp)*scale_tmp_zoo1;
        g_lim_pzres_sp = max(0,(1-exp(iv_pzres_sp*(thresh_pzres_sp-sp_v2))).*exp(-inh_szlzlp_sp*(lp_v2+sz_v2+lzres_v2)));
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_pzres_sp*t_lim_g_pzres_sp.*g_lim_pzres_sp.*oxy_lim_met.*pzres_v2)));
        
        % pzres grazing on lp
        t_lim_g_pzres_lp = exp(tc_g_pzres_lp*tmp)*scale_tmp_zoo1;
        g_lim_pzres_lp = max(0,(1-exp(iv_pzres_lp*(thresh_pzres_lp-lp_v2))).*exp(-inh_szlz_lp*(sz_v2+lzres_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_pzres_lp*t_lim_g_pzres_lp.*g_lim_pzres_lp.*oxy_lim_met.*pzres_v2)));
        
        % pzres Grazing on sz
        t_lim_g_pzres_sz = exp(tc_g_pzres_sz*tmp)*scale_tmp_zoo1;
        g_lim_pzres_sz = max(0,(1-exp(iv_pzres_sz*(thresh_pzres_sz-sz_v2))).*exp(-inh_lz_sz*lzres_v2));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_pzres_sz*t_lim_g_pzres_sz.*g_lim_pzres_sz.*oxy_lim_met.*pzres_v2)));
        
        % pzres Grazing on lzres
        t_lim_g_pzres_lz = exp(tc_g_pzres_lz*tmp)*scale_tmp_zoo1;
        g_lim_pzres_lz = max(0,1-exp(iv_pzres_lz*(thresh_pzres_lz-lzres_v2)));
        delta_lzres = lzres_v2-(lzres_v2./(1+dt./lzres_v2.*(gmax_pzres_lz*t_lim_g_pzres_lz.*g_lim_pzres_lz.*oxy_lim_met.*pzres_v2)));
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzres_v2 = lzres_v2-delta_lzres;
        pzres_v2 = pzres_v1+(ae_pzres-act_res_pzres)*(delta_sp+delta_lp+delta_sz+delta_lzres);
        NH_v2 = NH_v2+(act_res_pzres)*(delta_sp+delta_lp+delta_sz+delta_lzres);
        LPON_v2 = LPON_v2+(1-ae_pzres)*(delta_sp+delta_lp+delta_sz+delta_lzres);
        LOP_v2 = LOP_v2+delta_lp*r_SI_N;
        
        % Diagnostics
        SP2PZres=delta_sp;
        LP2PZres=delta_lp;
        SZ2PZres=delta_sz;
        SP2PZdvm=0;
        LP2PZdvm=0;
        SZ2PZdvm=0;
        LZres2PZres=delta_lzres;
        PZresgraz=delta_sp+delta_lp;
        PZdvmgraz=zeros(size(pzdvm_v2));
        LZdvm2PZres=0;
        LZres2PZdvm=0;
        LZdvm2PZdvm=0;
        PZgrazchl = delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;
        % -----------------------------%
        
        % *** END of Grazing Terms *** %
        
        
    else   %Night time grazing
        
        %----------------------------------------------%
        % --- Large Zooplankton (Resident) Grazing --- %
        %----------------------------------------------%
        
        % lzres grazing on sp
        t_lim_g_lzres_sp =exp(tc_g_lzres_sp*tmp)*scale_tmp_zoo1;
        g_lim_lzres_sp = max(0,1-exp(iv_lzres_sp*(thresh_lzres_sp-sp_v2)));
        oxy_lim_met = oxy_v2./(koxy_met+oxy_v2);
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_lzres_sp*t_lim_g_lzres_sp.*g_lim_lzres_sp.*oxy_lim_met.*lzres_v2)));
        
        % lzres grazing on lp
        t_lim_g_lzres_lp =exp(tc_g_lzres_lp*tmp)*scale_tmp_zoo1;
        g_lim_lzres_lp = max(0,1-exp(iv_lzres_lp*(thresh_lzres_lp-lp_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_lzres_lp*t_lim_g_lzres_lp.*g_lim_lzres_lp.*oxy_lim_met.*lzres_v2)));
        
        % lzres grazing on sz
        t_lim_g_lzres_sz = exp(tc_g_lzres_sz*tmp)*scale_tmp_zoo1;
        g_lim_lzres_sz = max(0,1-exp(iv_lzres_sz*(thresh_lzres_sz-sz_v2)));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_lzres_sz*t_lim_g_lzres_sz.*g_lim_lzres_sz.*oxy_lim_met.*lzres_v2)));
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzres_v2 = lzres_v1+(ae_lzres-act_res_lzres)*(delta_sp+delta_lp+delta_sz);
        NH_v2 = NH_v2+(act_res_lzres)*(delta_sp+delta_lp+delta_sz);
        LPON_v2 = LPON_v1+(1-ae_lzres)*(delta_sp+delta_lp+delta_sz);
        LOP_v2 = LOP_v1+delta_lp*r_SI_N;
        
        % Diagnostics
        SP2LZres=delta_sp;
        LP2LZres=delta_lp;
        SZ2LZres=delta_sz;
        LZresgraz = delta_sp+delta_lp;
        LZgrazchl = delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;  %Further modified below
        % -------------------------%
        
        
        %----------------------------------------------------------%
        % --- Large Zooplankton (Vertically-migrating) Grazing --- %
        %----------------------------------------------------------%
        
        % lzdvm grazing on sp
        t_lim_g_lzdvm_sp =exp(tc_g_lzdvm_sp*tmp)*scale_tmp_zoo1;
        g_lim_lzdvm_sp = max(0,1-exp(iv_lzdvm_sp*(thresh_lzdvm_sp-sp_v2)));
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_lzdvm_sp*t_lim_g_lzdvm_sp.*g_lim_lzdvm_sp.*oxy_lim_met.*lzdvm_v2)));
        
        % lzdvm grazing on lp
        t_lim_g_lzdvm_lp =exp(tc_g_lzdvm_lp*tmp)*scale_tmp_zoo1;
        g_lim_lzdvm_lp = max(0,1-exp(iv_lzdvm_lp*(thresh_lzdvm_lp-lp_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_lzdvm_lp*t_lim_g_lzdvm_lp.*g_lim_lzdvm_lp.*oxy_lim_met.*lzdvm_v2)));
        
        % lzdvm grazing on sz
        t_lim_g_lzdvm_sz = exp(tc_g_lzdvm_sz*tmp)*scale_tmp_zoo1;
        g_lim_lzdvm_sz = max(0,1-exp(iv_lzdvm_sz*(thresh_lzdvm_sz-sz_v2)));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_lzdvm_sz*t_lim_g_lzdvm_sz.*g_lim_lzdvm_sz.*oxy_lim_met.*lzdvm_v2)));
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzdvm_v2 = lzdvm_v1+(ae_lzdvm-act_res_lzdvm)*(delta_sp+delta_lp+delta_sz);
        NH_v2 = NH_v2+(act_res_lzdvm)*(delta_sp+delta_lp+delta_sz);
        LPON_v2 = LPON_v2+(1-ae_lzdvm)*(delta_sp+delta_lp+delta_sz);
        LOP_v2 = LOP_v2+delta_lp*r_SI_N;
        
        % Diagnostics
        SP2LZdvm=delta_sp;
        LP2LZdvm=delta_lp;
        SZ2LZdvm=delta_sz;
        LZdvmgraz = delta_sp+delta_lp;
        LZgrazchl = LZgrazchl + delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;
        % -------------------------%
        
        %--------------------------------------------------%
        % --- Predatory Zooplankton (Resident) Grazing --- %
        %--------------------------------------------------%
        
        % pzres grazing on sp
        lztot_v2 = lzres_v2 + lzdvm_v2;
        t_lim_g_pzres_sp = exp(tc_g_pzres_sp*tmp)*scale_tmp_zoo1;
        g_lim_pzres_sp = max(0,(1-exp(iv_pzres_sp*(thresh_pzres_sp-sp_v2))).*exp(-inh_szlzlp_sp*(lp_v2+sz_v2+lztot_v2)));
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_pzres_sp*t_lim_g_pzres_sp.*g_lim_pzres_sp.*oxy_lim_met.*pzres_v2)));
        
        % pzres grazing on lp
        t_lim_g_pzres_lp = exp(tc_g_pzres_lp*tmp)*scale_tmp_zoo1;
        g_lim_pzres_lp = max(0,(1-exp(iv_pzres_lp*(thresh_pzres_lp-lp_v2))).*exp(-inh_szlz_lp*(sz_v2+lztot_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_pzres_lp*t_lim_g_pzres_lp.*g_lim_pzres_lp.*oxy_lim_met.*pzres_v2)));
        
        % pzres Grazing on sz
        t_lim_g_pzres_sz = exp(tc_g_pzres_sz*tmp)*scale_tmp_zoo1;
        g_lim_pzres_sz = max(0,(1-exp(iv_pzres_sz*(thresh_pzres_sz-sz_v2))).*exp(-inh_lz_sz*lztot_v2));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_pzres_sz*t_lim_g_pzres_sz.*g_lim_pzres_sz.*oxy_lim_met.*pzres_v2)));
        
        % pzres Grazing on lzres and lzdvm
        t_lim_g_pzres_lz = exp(tc_g_pzres_lz*tmp)*scale_tmp_zoo1;
        g_lim_pzres_lz = max(0,1-exp(iv_pzres_lz*(thresh_pzres_lz-lztot_v2)));
        delta_lztot = lztot_v2-(lztot_v2./(1+dt./lztot_v2.*(gmax_pzres_lz*t_lim_g_pzres_lz.*g_lim_pzres_lz.*oxy_lim_met.*pzres_v2)));
        delta_lzres = delta_lztot .* lzres_v2./lztot_v2;
        delta_lzdvm = delta_lztot .* lzdvm_v2./lztot_v2;
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzres_v2 = lzres_v2-delta_lzres;
        lzdvm_v2 = lzdvm_v2-delta_lzdvm;
        pzres_v2 = pzres_v1+(ae_pzres-act_res_pzres)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        NH_v2 = NH_v2+(act_res_pzres)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        LPON_v2 = LPON_v2+(1-ae_pzres)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        LOP_v2 = LOP_v2+delta_lp*r_SI_N;
        
        % Diagnostics
        SP2PZres=delta_sp;
        LP2PZres=delta_lp;
        SZ2PZres=delta_sz;
        LZres2PZres=delta_lzres;
        LZdvm2PZres=delta_lzdvm;
        PZresgraz=delta_sp+delta_lp;
        PZgrazchl = delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;   %Further modified below
        
        %---------------------------------------------------%
        % --- Predatory Zooplankton (Migratory) Grazing --- %
        %---------------------------------------------------%
        
        % pzdvm grazing on sp
        lztot_v2 = lzres_v2 + lzdvm_v2;
        t_lim_g_pzdvm_sp = exp(tc_g_pzdvm_sp*tmp)*scale_tmp_zoo1;
        g_lim_pzdvm_sp = max(0,(1-exp(iv_pzdvm_sp*(thresh_pzdvm_sp-sp_v2))).*exp(-inh_szlzlp_sp*(lp_v2+sz_v2+lztot_v2)));
        delta_sp = sp_v2-(sp_v2./(1+dt./sp_v2.*(gmax_pzdvm_sp*t_lim_g_pzdvm_sp.*g_lim_pzdvm_sp.*oxy_lim_met.*pzdvm_v2)));
        
        % pzdvm grazing on lp
        t_lim_g_pzdvm_lp = exp(tc_g_pzdvm_lp*tmp)*scale_tmp_zoo1;
        g_lim_pzdvm_lp = max(0,(1-exp(iv_pzdvm_lp*(thresh_pzdvm_lp-lp_v2))).*exp(-inh_szlz_lp*(sz_v2+lztot_v2)));
        delta_lp = lp_v2-(lp_v2./(1+dt./lp_v2.*(gmax_pzdvm_lp*t_lim_g_pzdvm_lp.*g_lim_pzdvm_lp.*oxy_lim_met.*pzdvm_v2)));
        
        % pzdvm Grazing on sz
        t_lim_g_pzdvm_sz = exp(tc_g_pzdvm_sz*tmp)*scale_tmp_zoo1;
        g_lim_pzdvm_sz = max(0,(1-exp(iv_pzdvm_sz*(thresh_pzdvm_sz-sz_v2))).*exp(-inh_lz_sz*lztot_v2));
        delta_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(gmax_pzdvm_sz*t_lim_g_pzdvm_sz.*g_lim_pzdvm_sz.*oxy_lim_met.*pzdvm_v2)));
        
        % pzdvm Grazing on lzres and lzdvm
        t_lim_g_pzdvm_lz = exp(tc_g_pzdvm_lz*tmp)*scale_tmp_zoo1;
        g_lim_pzdvm_lz = max(0,1-exp(iv_pzdvm_lz*(thresh_pzdvm_lz-lztot_v2)));
        delta_lztot = lztot_v2-(lztot_v2./(1+dt./lztot_v2.*(gmax_pzdvm_lz*t_lim_g_pzdvm_lz.*g_lim_pzdvm_lz.*oxy_lim_met.*pzdvm_v2)));
        delta_lzres = delta_lztot .* lzres_v2./lztot_v2;
        delta_lzdvm = delta_lztot .* lzdvm_v2./lztot_v2;
        
        % Update
        sp_v2 = sp_v2-delta_sp;
        lp_v2 = lp_v2-delta_lp;
        sz_v2 = sz_v2-delta_sz;
        lzres_v2 = lzres_v2-delta_lzres;
        lzdvm_v2 = lzdvm_v2-delta_lzdvm;
        pzdvm_v2 = pzdvm_v1+(ae_pzdvm-act_res_pzdvm)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        NH_v2 = NH_v2+(act_res_pzdvm)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        LPON_v2 = LPON_v2+(1-ae_pzdvm)*(delta_sp+delta_lp+delta_sz+delta_lztot);
        LOP_v2 = LOP_v2+delta_lp*r_SI_N;
        
        % Diagnostics
        SP2PZdvm=delta_sp;
        LP2PZdvm=delta_lp;
        SZ2PZdvm=delta_sz;
        LZres2PZdvm=delta_lzres;
        LZdvm2PZdvm=delta_lzdvm;
        PZdvmgraz=delta_sp+delta_lp;
        PZgrazchl = PZgrazchl + delta_sp * carbon_nitrogen * 12 .* chl2c_sp + delta_lp * carbon_nitrogen * 12 .* chl2c_lp;
        
    end   %End if statement for day/night grazing
    
    
    
    if day==1 & euponly==1
        
        %-------------------------------------%
        % *** Zooplankton Mortality Terms *** %
        %-------------------------------------%
        
        % sz Mortality
        t_lim_m_sz = exp(tc_m_sz*tmp)*scale_tmp_zoo2;
        mortality_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(ref_mort_sz*t_lim_m_sz.*sz_v2)));
        
        % lzres Mortality
        t_lim_m_lzres = exp(tc_m_lz*tmp)*scale_tmp_zoo2;
        mortality_lzres = lzres_v2-(lzres_v2./(1+dt./lzres_v2.*(mort_day_lzres*t_lim_m_lzres.*oxy_lim_met.*lzres_v2.^2)));
        
        % lzdvm Mortality
        t_lim_m_lzdvm = exp(tc_m_lz*tmpdeep)*scale_tmp_zoo2;
        mortality_lzdvm = lzdvm_v2-(lzdvm_v2./(1+dt./lzdvm_v2.*(mort_day_lzdvm*t_lim_m_lzdvm.*oxy_lim_met.*lzdvm_v2.^2)));
        
        % pzdvm Mortality
        t_lim_m_pzdvm = exp(tc_m_pz*tmpdeep)*scale_tmp_zoo2;
        mortality_pzdvm = pzdvm_v2-(pzdvm_v2./(1+dt./pzdvm_v2.*(mort_day_pzdvm*t_lim_m_pzdvm.*oxy_lim_met.*pzdvm_v2.^2)));
        
        % pzdvm Mortality
        t_lim_m_pzres = exp(tc_m_pz*tmp)*scale_tmp_zoo2;
        mortality_pzres = pzres_v2-(pzres_v2./(1+dt./pzres_v2.*(mort_day_pzres*t_lim_m_pzres.*oxy_lim_met.*pzres_v2.^2)));
        
        % Update
        sz_v2 = sz_v2-mortality_sz;
        lzres_v2 = lzres_v2-mortality_lzres;
        pzres_v2 = pzres_v2-mortality_pzres;
        lzdvm_v2 = lzdvm_v1-mortality_lzdvm;
        pzdvm_v2 = pzdvm_v1-mortality_pzdvm;
        LPON_v2 = LPON_v2 + mortality_sz + mortality_lzres + mortality_pzres;
        %LPON_v2 = LPON_v2 + mortality_lzdvm + mortality_pzdvm;           %---------only turn on to test conservation
        DVM_mortality=mortality_lzdvm + mortality_pzdvm;
        % ----------------- %
        
        
        %-----------------------------------------------%
        % ******* Zooplankton Excretion Terms ********* %
        % Excretion = exp(a0 + a2 * ln(NW)) * Temp^a2
        %-----------------------------------------------%
        
        % lzres Excretion
        t_lim_m_lzres = exp(Ikeda_a2*tmp);
        basalexcretion_lzres = lzres_v2-(lzres_v2./(1+dt./lzres_v2.*(Ikeda_lz*t_lim_m_lzres.*oxy_lim_met.*lzres_v2)));
        
        % lzdvm Excretion
        t_lim_m_lzdvm = exp(Ikeda_a2*tmpdeep);
        basalexcretion_lzdvm = lzdvm_v2-(lzdvm_v2./(1+dt./lzdvm_v2.*(Ikeda_lz*t_lim_m_lzdvm.*oxy_lim_met.*lzdvm_v2)));
        
        % pzres Excretion
        t_lim_m_pzres = exp(Ikeda_a2*tmp);
        basalexcretion_pzres = pzres_v2-(pzres_v2./(1+dt./pzres_v2.*(Ikeda_pz*t_lim_m_pzres.*oxy_lim_met.*pzres_v2)));
        
        % pzdvm Excretion
        t_lim_m_pzdvm = exp(Ikeda_a2*tmpdeep);
        basalexcretion_pzdvm = pzdvm_v2-(pzdvm_v2./(1+dt./pzdvm_v2.*(Ikeda_pz*t_lim_m_pzdvm.*oxy_lim_met.*pzdvm_v2)));
        
        % Update
        lzres_v2 = lzres_v2-basalexcretion_lzres;
        pzres_v2 = pzres_v2-basalexcretion_pzres;
        lzdvm_v2 = lzdvm_v2-basalexcretion_lzdvm;
        pzdvm_v2 = pzdvm_v2-basalexcretion_pzdvm;
        NH_v2 = NH_v2 + basalexcretion_lzres + basalexcretion_pzres;
        %NH_v2 = NH_v2 + basalexcretion_lzdvm + basalexcretion_pzdvm;           %---------only turn on to test conservation
        DVM_excretion=basalexcretion_lzdvm + basalexcretion_pzdvm;
        % ----------------- %
        
    else    %Nighttime
        
        %-------------------------------------%
        % *** Zooplankton Mortality Terms *** %
        %-------------------------------------%
        
        if day==1
            %The original code was written for a euphotic zone only model - this code section allows for a model that explicitly includes mesozooplankton vertical migration between the euphotic zone and the mesopelagic, in which case mesozooplankton need to undergo all of their activities during the day (as well as at night), but with rates that are appropriate for daytime activity.
            mort_night_lzres=mort_day_lzres;
            mort_night_lzdvm=mort_day_lzdvm;
            mort_night_pzres=mort_day_pzres;
            mort_night_pzdvm=mort_day_pzdvm;
        end
        
        % sz Mortality
        t_lim_m_sz = exp(tc_m_sz*tmp)*scale_tmp_zoo2;
        mortality_sz = sz_v2-(sz_v2./(1+dt./sz_v2.*(ref_mort_sz*t_lim_m_sz.*sz_v2)));
        
        % lzres Mortality
        t_lim_m_lzres = exp(tc_m_lz*tmp)*scale_tmp_zoo2;
        mortality_lzres = lzres_v2-(lzres_v2./(1+dt./lzres_v2.*(mort_night_lzres*t_lim_m_lzres.*oxy_lim_met.*lzres_v2.^2)));
        
        % lzdvm Mortality
        t_lim_m_lzdvm = exp(tc_m_lz*tmp)*scale_tmp_zoo2;
        mortality_lzdvm = lzdvm_v2-(lzdvm_v2./(1+dt./lzdvm_v2.*(mort_night_lzdvm*t_lim_m_lzdvm.*oxy_lim_met.*lzdvm_v2.^2)));
        
        % pzres Mortality
        t_lim_m_pzres = exp(tc_m_pz*tmp)*scale_tmp_zoo2;
        mortality_pzres = pzres_v2-(pzres_v2./(1+dt./pzres_v2.*(mort_night_pzres*t_lim_m_pzres.*oxy_lim_met.*pzres_v2.^2)));
        
        % pzdvm Mortality
        t_lim_m_pzdvm = exp(tc_m_pz*tmp)*scale_tmp_zoo2;
        mortality_pzdvm = pzdvm_v2-(pzdvm_v2./(1+dt./pzdvm_v2.*(mort_night_pzdvm*t_lim_m_pzdvm.*oxy_lim_met.*pzdvm_v2.^2)));
        
        % Update
        sz_v2 = sz_v2-mortality_sz;
        lzres_v2 = lzres_v2-mortality_lzres;
        pzres_v2 = pzres_v2-mortality_pzres;
        lzdvm_v2 = lzdvm_v2-mortality_lzdvm;
        pzdvm_v2 = pzdvm_v2-mortality_pzdvm;
        LPON_v2 = LPON_v2 + mortality_sz + mortality_lzres + mortality_pzres + mortality_lzdvm + mortality_pzdvm;
        DVM_excretion=zeros(size(lzres_v2));
        % ----------------- %
        
        %-------------------------------------%
        % *** Zooplankton Excretion Terms *** %
        %-------------------------------------%
        
        % lzres Excretion
        t_lim_m_lzres = exp(Ikeda_a2*tmp);
        basalexcretion_lzres = lzres_v2-(lzres_v2./(1+dt./lzres_v2.*(Ikeda_lz*t_lim_m_lzres.*oxy_lim_met.*lzres_v2)));
        
        % lzdvm Excretion
        t_lim_m_lzdvm = exp(Ikeda_a2*tmp);
        basalexcretion_lzdvm = lzdvm_v2-(lzdvm_v2./(1+dt./lzdvm_v2.*(Ikeda_lz*t_lim_m_lzdvm.*oxy_lim_met.*lzdvm_v2)));
        
        % pzres Excretion
        t_lim_m_pzres = exp(Ikeda_a2*tmp);
        basalexcretion_pzres = pzres_v2-(pzres_v2./(1+dt./pzres_v2.*(Ikeda_pz*t_lim_m_pzres.*oxy_lim_met.*pzres_v2)));
        
        % pzdvm Excretion
        t_lim_m_pzdvm = exp(Ikeda_a2*tmp);
        basalexcretion_pzdvm = pzdvm_v2-(pzdvm_v2./(1+dt./pzdvm_v2.*(Ikeda_pz*t_lim_m_pzdvm.*oxy_lim_met.*pzdvm_v2)));
        
        % Update
        lzres_v2 = lzres_v2-basalexcretion_lzres;
        pzres_v2 = pzres_v2-basalexcretion_pzres;
        lzdvm_v2 = lzdvm_v2-basalexcretion_lzdvm;
        pzdvm_v2 = pzdvm_v2-basalexcretion_pzdvm;
        NH_v2 = NH_v2 + basalexcretion_lzres + basalexcretion_pzres + basalexcretion_lzdvm + basalexcretion_pzdvm;
        DVM_mortality=zeros(size(lzres_v2));
        % ----------------- %
        % *** END of Mortality Terms *** %
        
    end
    
    % ### END Functional Groups ### %
    
    %----------------------------------%
    % ### Chemical Transformations ### %
    %----------------------------------%
    
    % Nitrification
    t_lim_nit =exp(tc_nitr*tmp)*scale_tmp_decomp;
    nitrification_NH = NH_v2-(NH_v2./(1+dt.*(ref_nitr*t_lim_nit.*oxy_lim_mic)));
    NH_v2 = NH_v2-nitrification_NH;
    NO_v2 = NO_v2+nitrification_NH;
    
    % Aggregation of PON to LPON
    t_lim_dec_PON_LPON = exp(tc_dec_PON_LPON*tmp)*scale_tmp_decomp;
    aggregation_PON = PON_v2-(PON_v2./(1+dt.*(ref_dec_PON_LPON*t_lim_dec_PON_LPON.*PON_v2)));
    PON_v2 = PON_v2-aggregation_PON;
    LPON_v2 = LPON_v2+aggregation_PON;
    
    % Decomp of PON to NH
    t_lim_dec_PON_NH = exp(tc_dec_PON_NH*tmp)*scale_tmp_decomp;
    remineralization_PON = PON_v2-(PON_v2./(1+dt.*(ref_dec_PON_NH*t_lim_dec_PON_NH.*oxy_lim_mic)));
    PON_v2 = PON_v2-remineralization_PON;
    NH_v2 = NH_v2+remineralization_PON;
    
    % Decomp of LPON to NH
    t_lim_dec_PON_NH = exp(tc_dec_PON_NH*tmp)*scale_tmp_decomp;
    remineralization_LPON = LPON_v2-(LPON_v2./(1+dt.*(ref_dec_LPON_NH*t_lim_dec_PON_NH.*oxy_lim_mic)));
    LPON_v2 = LPON_v2-remineralization_LPON;
    NH_v2 = NH_v2+remineralization_LPON;
    
    % Decomp of DON to NH
    t_lim_dec_DON_NH = exp(tc_dec_DON_NH*tmp)*scale_tmp_decomp;
    remineralization_DON = DON_v2-(DON_v2./(1+dt.*(ref_dec_DON_NH*t_lim_dec_DON_NH.*oxy_lim_mic)));
    DON_v2 = DON_v2-remineralization_DON;
    NH_v2 = NH_v2+remineralization_DON*(1-f_DON_DONref);
    DONref_v2 = DONref_v1+remineralization_DON*f_DON_DONref;
        
    % Decomp of Refractory DON to NH
    t_lim_dec_DON_NH = exp(tc_dec_DON_NH*tmp)*scale_tmp_decomp;
    remineralization_DONref = DONref_v2-(DONref_v2./(1+dt.*(ref_dec_DONref_NH*t_lim_dec_DON_NH.*oxy_lim_mic)));
    DONref_v2 = DONref_v2-remineralization_DONref;
    NH_v2 = NH_v2+remineralization_DONref;
    
    % Decomp of PON to DON
    t_lim_dec_PON_DON = exp(tc_dec_PON_DON*tmp)*scale_tmp_decomp;
    dissolution_PON = PON_v2-(PON_v2./(1+dt.*(ref_dec_PON_DON*t_lim_dec_PON_DON.*oxy_lim_mic)));
    PON_v2 = PON_v2-dissolution_PON;
    DON_v2 = DON_v2+dissolution_PON;
    
    % Decomp of LPON to DON
    t_lim_dec_PON_DON = exp(tc_dec_PON_DON*tmp)*scale_tmp_decomp;
    dissolution_LPON = LPON_v2-(LPON_v2./(1+dt.*(ref_dec_LPON_DON*t_lim_dec_PON_DON.*oxy_lim_mic)));
    LPON_v2 = LPON_v2-dissolution_LPON;
    DON_v2 = DON_v2+dissolution_LPON;
    
    % Decomp of OP to SI (First OP Loss Term)
    t_lim_dec_OP_SI = exp(tc_dec_OP_SI*tmp)*scale_tmp_decomp;
    delta_OP = OP_v2-(OP_v2./(1+dt.*(ref_dec_OP_SI*t_lim_dec_OP_SI)));
    delta_OP = max(0,delta_OP);
    OP_v2 = OP_v2-delta_OP;
    SI_v2 = SI_v2+delta_OP;
    
    % Decomp of LOP to SI (First OP Loss Term)
    t_lim_dec_OP_SI = exp(tc_dec_OP_SI*tmp)*scale_tmp_decomp;
    delta_LOP = LOP_v2-(LOP_v2./(1+dt.*(ref_dec_OP_SI*t_lim_dec_OP_SI)));
    LOP_old = LOP_v2;
    LOP_v2 = LOP_v2-delta_LOP;
    SI_v2 = SI_v2+delta_LOP;
    
    % ### END Chemical Transformations ### %
    
end % Iter for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---------------------------------------------------------------------%
%-------------------Update Oxygen and DIC-----------------------------%
%---------------------------------------------------------------------%
if CO2on==1
    OrgNBegin = sp_v1 + lp_v1 + sz_v1 + lzres_v1 + lzdvm_v1 + pzres_v1 + ...
        pzdvm_v1 + PON_v1 + LPON_v1 + DON_v1 + DONref_v1;
    OrgNEnd = sp_v2 + lp_v2 + sz_v2 + lzres_v2 + lzdvm_v2 + pzres_v2 + ...
        pzdvm_v2 + PON_v2 + LPON_v2 + DON_v2 + DONref_v2;
    OrgNChange = OrgNEnd - OrgNBegin;
    diff = OrgNChange * oxy_nitrogen_nh4 ...
        - (NO_v2 - NO_v1) * (oxy_nitrogen_no3 - oxy_nitrogen_nh4); 
    for j=1:length(diff)
        if diff(j)>0
            oxy_v2(j) = oxy_v1(j) + diff(j) ;
        else
            oxy_v2(j) = oxy_v1(j) ./ (1-diff(j)./oxy_v2(j));  %A backwards Euler implementation for oxygen to be used when O2 is decreasing to avoid negative overshoots
        end
    end
    DIC_v2 = DIC_v1 - OrgNChange * carbon_nitrogen;
    alk_v2 = alk_v1 - (NO_v2 - NO_v1) * carbon_nitrogen + (NH_v2 - NH_v1)  ...
        - (NO_v2 - NO_v1 + NH_v2 - NH_v1)  / nitrogen_phosphorus;
end


if Thoriumon==1
    kforward = 0.013*6.625;  %thorium adsorption constant = 0.13 m3 / mmol C / day, from Stukel et al. (2019)
    kbackward = 3/365.25;    %thorium desorption constant = 3 yr-1, from Resplandy et al. (2012)
    LZfactor = 0.04;         %This is a factor that accounts for lower expected adsorption onto mesozooplankton (relative to protists and detritus) as a result of a lower surface area:volume ratio
    PZfactor = 0.01;         %This is a factor that accounts for lower expected adsorption onto mesozooplankton (relative to protists and detritus) as a result of a lower surface area:volume ratio
    lam234 = 0.0287612937991679;  %Th-234 decay constant (d-1)
    
    dTh_v1=tracer(:,iThorium+0);        dTh_v2=dTh_v1;
    spTh_v1=tracer(:,iThorium+1);       spTh_v2=spTh_v1;
    lpTh_v1=tracer(:,iThorium+2);       lpTh_v2=lpTh_v1;
    szTh_v1=tracer(:,iThorium+3);       szTh_v2=szTh_v1;
    lzresTh_v1=tracer(:,iThorium+4);    lzresTh_v2=lzresTh_v1;
    lzdvmTh_v1=tracer(:,iThorium+5);    lzdvmTh_v2=lzdvmTh_v1;
    pzresTh_v1=tracer(:,iThorium+6);    pzresTh_v2=pzresTh_v1;
    pzdvmTh_v1=tracer(:,iThorium+7);    pzdvmTh_v2=pzdvmTh_v1;
    ponTh_v1=tracer(:,iThorium+8);      ponTh_v2=ponTh_v1;
    lponTh_v1=tracer(:,iThorium+9);     lponTh_v2=lponTh_v1;
    
    %------------Grazing thorium transformations-----------------------------------------------------------------------------------------
    spTh_v2 = spTh_v2 - SP2SZ.*spTh_v1./sp_v1;
    szTh_v2 = szTh_v2 + SP2SZ.*gge_sz.*spTh_v1./sp_v1;
    ponTh_v2 = ponTh_v2 + SP2SZ.*(1-ae_sz).*spTh_v1./sp_v1;
    dTh_v2 = dTh_v2 + SP2SZ.*(ae_sz-gge_sz).*spTh_v1./sp_v1;
    
    lpTh_v2 = lpTh_v2 - LP2SZ.*lpTh_v1./lp_v1;
    szTh_v2 = szTh_v2 + LP2SZ.*gge_sz.*lpTh_v1./lp_v1;
    ponTh_v2 = ponTh_v2 + LP2SZ.*(1-ae_sz).*lpTh_v1./lp_v1;
    dTh_v2 = dTh_v2 + LP2SZ.*(ae_sz-gge_sz).*lpTh_v1./lp_v1;
    
    spTh_v2 = spTh_v2 - SP2LZres.*spTh_v1./sp_v1;
    lponTh_v2 = lponTh_v2 + SP2LZres.*(1-act_res_lzres).*spTh_v1./sp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SP2LZres.*act_res_lzres.*spTh_v1./sp_v1;
    
    spTh_v2 = spTh_v2 - SP2LZdvm.*spTh_v1./sp_v1;
    lponTh_v2 = lponTh_v2 + SP2LZdvm.*(1-act_res_lzdvm).*spTh_v1./sp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SP2LZdvm.*act_res_lzdvm.*spTh_v1./sp_v1;
    
    spTh_v2 = spTh_v2 - SP2PZres.*spTh_v1./sp_v1;
    lponTh_v2 = lponTh_v2 + SP2PZres.*(1-act_res_pzres).*spTh_v1./sp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SP2PZres.*act_res_pzres.*spTh_v1./sp_v1;
    
    spTh_v2 = spTh_v2 - SP2PZdvm.*spTh_v1./sp_v1;
    lponTh_v2 = lponTh_v2 + SP2PZdvm.*(1-act_res_pzdvm).*spTh_v1./sp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SP2PZdvm.*act_res_pzdvm.*spTh_v1./sp_v1;
    
    lpTh_v2 = lpTh_v2 - LP2LZres.*lpTh_v1./lp_v1;
    lponTh_v2 = lponTh_v2 + LP2LZres.*(1-act_res_lzres).*lpTh_v1./lp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LP2LZres.*act_res_lzres.*lpTh_v1./lp_v1;
    
    lpTh_v2 = lpTh_v2 - LP2LZdvm.*lpTh_v1./lp_v1;
    lponTh_v2 = lponTh_v2 + LP2LZdvm.*(1-act_res_lzdvm).*lpTh_v1./lp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LP2LZdvm.*act_res_lzdvm.*lpTh_v1./lp_v1;
    
    lpTh_v2 = lpTh_v2 - LP2PZres.*lpTh_v1./lp_v1;
    lponTh_v2 = lponTh_v2 + LP2PZres.*(1-act_res_pzres).*lpTh_v1./lp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LP2PZres.*act_res_pzres.*lpTh_v1./lp_v1;

    lpTh_v2 = lpTh_v2 - LP2PZdvm.*lpTh_v1./lp_v1;
    lponTh_v2 = lponTh_v2 + LP2PZdvm.*(1-act_res_pzdvm).*lpTh_v1./lp_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LP2PZdvm.*act_res_pzdvm.*lpTh_v1./lp_v1;

    szTh_v2 = szTh_v2 - SZ2LZres.*szTh_v1./sz_v1;
    lponTh_v2 = lponTh_v2 + SZ2LZres.*(1-act_res_lzres).*szTh_v1./sz_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SZ2LZres.*act_res_lzres.*szTh_v1./sz_v1;
    
    szTh_v2 = szTh_v2 - SZ2LZdvm.*szTh_v1./sz_v1;
    lponTh_v2 = lponTh_v2 + SZ2LZdvm.*(1-act_res_lzdvm).*szTh_v1./sz_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SZ2LZdvm.*act_res_lzdvm.*szTh_v1./sz_v1;
    
    szTh_v2 = szTh_v2 - SZ2PZres.*szTh_v1./sz_v1;
    lponTh_v2 = lponTh_v2 + SZ2PZres.*(1-act_res_pzres).*szTh_v1./sz_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SZ2PZres.*act_res_pzres.*szTh_v1./sz_v1;
    
    szTh_v2 = szTh_v2 - SZ2PZdvm.*szTh_v1./sz_v1;
    lponTh_v2 = lponTh_v2 + SZ2PZdvm.*(1-act_res_pzdvm).*szTh_v1./sz_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + SZ2PZdvm.*act_res_pzdvm.*szTh_v1./sz_v1;
    
    lzresTh_v2 = lzresTh_v2 - LZres2PZres.*lzresTh_v1./lzres_v1;
    lponTh_v2 = lponTh_v2 + LZres2PZres.*(1-act_res_pzres).*lzresTh_v1./lzres_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LZres2PZres.*act_res_pzres.*lzresTh_v1./lzres_v1;
    
    lzdvmTh_v2 = lzdvmTh_v2 - LZdvm2PZres.*lzdvmTh_v1./lzdvm_v1;
    lponTh_v2 = lponTh_v2 + LZdvm2PZres.*(1-act_res_pzres).*lzdvmTh_v1./lzdvm_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LZdvm2PZres.*act_res_pzres.*lzdvmTh_v1./lzdvm_v1;
    
    lzresTh_v2 = lzresTh_v2 - LZres2PZdvm.*lzresTh_v1./lzres_v1;
    lponTh_v2 = lponTh_v2 + LZres2PZdvm.*(1-act_res_pzdvm).*lzresTh_v1./lzres_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LZres2PZdvm.*act_res_pzdvm.*lzresTh_v1./lzres_v1;
    
    lzdvmTh_v2 = lzdvmTh_v2 - LZdvm2PZdvm.*lzdvmTh_v1./lzdvm_v1;
    lponTh_v2 = lponTh_v2 + LZdvm2PZdvm.*(1-act_res_pzdvm).*lzdvmTh_v1./lzdvm_v1;   %Note, there is evidence that fecal pellets are enriched in thorium relative to zoop prey, and zoops have low thorium, so I assume that thorium that would go to biomass goes to fecal pellets instead
    dTh_v2 = dTh_v2 + LZdvm2PZdvm.*act_res_pzdvm.*lzdvmTh_v1./lzdvm_v1;
    
    
    %------------Mortality thorium transformations-----------------------------------------------------------------------------------------
    spTh_v2 = spTh_v2 - mortality_sp.*spTh_v1./sp_v1;
    ponTh_v2 = ponTh_v2 + mortality_sp.*spTh_v1./sp_v1;
    
    lpTh_v2 = lpTh_v2 - mortality_lp.*lpTh_v1./lp_v1;
    ponTh_v2 = ponTh_v2 + mortality_lp.*lpTh_v1./lp_v1;
    
    szTh_v2 = szTh_v2 - mortality_sz.*szTh_v1./sz_v1;
    lponTh_v2 = lponTh_v2 + mortality_sz.*szTh_v1./sz_v1;
    
    lzresTh_v2 = lzresTh_v2 - mortality_lzres.*lzresTh_v1./lzres_v1;
    lponTh_v2 = lponTh_v2 + mortality_lzres.*lzresTh_v1./lzres_v1;
    
    lzdvmTh_v2 = lzdvmTh_v2 - mortality_lzdvm.*lzdvmTh_v1./lzdvm_v1;
    lponTh_v2 = lponTh_v2 + mortality_lzdvm.*lzdvmTh_v1./lzdvm_v1;
    
    pzresTh_v2 = pzresTh_v2 - mortality_pzres.*pzresTh_v1./pzres_v1;
    lponTh_v2 = lponTh_v2 + mortality_pzres.*pzresTh_v1./pzres_v1;
    
    pzdvmTh_v2 = pzdvmTh_v2 - mortality_pzdvm.*pzdvmTh_v1./pzdvm_v1;
    lponTh_v2 = lponTh_v2 + mortality_pzdvm.*pzdvmTh_v1./pzdvm_v1;
    
    %------------Respiration and Excretion-----------------------------------------------------------------------------------------
    spTh_v2 = spTh_v2 - respiration_sp.*spTh_v1./sp_v1;
    dTh_v2 = dTh_v2 + respiration_sp.*spTh_v1./sp_v1;
    
    lpTh_v2 = lpTh_v2 - respiration_lp.*lpTh_v1./lp_v1;
    dTh_v2 = dTh_v2 + respiration_lp.*lpTh_v1./lp_v1;
    
    lzresTh_v2 = lzresTh_v2 - basalexcretion_lzres.*lzresTh_v1./lzres_v1;
    dTh_v2 = dTh_v2 + basalexcretion_lzres.*lzresTh_v1./lzres_v1;
    
    lzdvmTh_v2 = lzdvmTh_v2 - basalexcretion_lzdvm.*lzdvmTh_v1./lzdvm_v1;
    dTh_v2 = dTh_v2 + basalexcretion_lzdvm.*lzdvmTh_v1./lzdvm_v1;
    
    pzresTh_v2 = pzresTh_v2 - basalexcretion_pzres.*pzresTh_v1./pzres_v1;
    dTh_v2 = dTh_v2 + basalexcretion_pzres.*pzresTh_v1./pzres_v1;
    
    pzdvmTh_v2 = pzdvmTh_v2 - basalexcretion_pzdvm.*pzdvmTh_v1./pzdvm_v1;
    dTh_v2 = dTh_v2 + basalexcretion_pzdvm.*pzdvmTh_v1./pzdvm_v1;
    
    %------------Aggregation-----------------------------------------------------------------------------------------
    ponTh_v2 = ponTh_v2 - aggregation_PON.*ponTh_v1./PON_v1;
    lponTh_v2 = lponTh_v2 + aggregation_PON.*ponTh_v1./PON_v1;
    
    %------------Remineralization-----------------------------------------------------------------------------------------
    ponTh_v2 = ponTh_v2 - remineralization_PON.*ponTh_v1./PON_v1;
    dTh_v2 = dTh_v2 + remineralization_PON.*ponTh_v1./PON_v1;
    
    lponTh_v2 = lponTh_v2 - remineralization_LPON.*lponTh_v1./LPON_v1;
    dTh_v2 = dTh_v2 + remineralization_LPON.*lponTh_v1./LPON_v1;
    
    
    ponTh_v2 = ponTh_v2 - dissolution_PON.*ponTh_v1./PON_v1;
    dTh_v2 = dTh_v2 + dissolution_PON.*ponTh_v1./PON_v1;
    
    lponTh_v2 = lponTh_v2 - dissolution_LPON.*lponTh_v1./LPON_v1;
    dTh_v2 = dTh_v2 + dissolution_LPON.*lponTh_v1./LPON_v1;
    
    %------------Desorption----------------------------------------------------------------------------------------------
    spTh_v2 = spTh_v2 - spTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + spTh_v1.*kbackward*dt;
    
    lpTh_v2 = lpTh_v2 - lpTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + lpTh_v1.*kbackward*dt;
    
    szTh_v2 = szTh_v2 - szTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + szTh_v1.*kbackward*dt;
    
    lzresTh_v2 = lzresTh_v2 - lzresTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + lzresTh_v1.*kbackward*dt;
    
    lzdvmTh_v2 = lzdvmTh_v2 - lzdvmTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + lzdvmTh_v1.*kbackward*dt;
    
    pzresTh_v2 = pzresTh_v2 - pzresTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + pzresTh_v1.*kbackward*dt;
    
    pzdvmTh_v2 = pzdvmTh_v2 - pzdvmTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + pzdvmTh_v1.*kbackward*dt;
    
    ponTh_v2 = ponTh_v2 - ponTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + ponTh_v1.*kbackward*dt;
    
    lponTh_v2 = lponTh_v2 - lponTh_v1.*kbackward*dt;
    dTh_v2 = dTh_v2 + lponTh_v1.*kbackward*dt;
    
    %------------Sorption----------------------------------------------------------------------------------------------
    spTh_v2 = spTh_v2 + sp_v1.*dTh_v1.*kforward*dt;
    dTh_v2 = dTh_v2 - sp_v1.*dTh_v1.*kforward*dt;
    
    lpTh_v2 = lpTh_v2 + lp_v1.*dTh_v1.*kforward*dt;
    dTh_v2 = dTh_v2 - lp_v1.*dTh_v1.*kforward*dt;
    
    szTh_v2 = szTh_v2 + sz_v1.*dTh_v1.*kforward*dt;
    dTh_v2 = dTh_v2 - sz_v1.*dTh_v1.*kforward*dt;
    
    lzresTh_v2 = lzresTh_v2 + lzres_v1.*dTh_v1.*kforward.*LZfactor*dt;
    dTh_v2 = dTh_v2 - lzres_v1.*dTh_v1.*kforward.*LZfactor*dt;
    
    lzdvmTh_v2 = lzdvmTh_v2 + lzdvm_v1.*dTh_v1.*kforward.*LZfactor*dt;
    dTh_v2 = dTh_v2 - lzdvm_v1.*dTh_v1.*kforward.*LZfactor*dt;
    
    pzresTh_v2 = pzresTh_v2 + pzres_v1.*dTh_v1.*kforward.*PZfactor*dt;
    dTh_v2 = dTh_v2 - pzres_v1.*dTh_v1.*kforward.*PZfactor*dt;
    
    pzdvmTh_v2 = pzdvmTh_v2 + pzdvm_v1.*dTh_v1.*kforward.*PZfactor*dt;
    dTh_v2 = dTh_v2 - pzdvm_v1.*dTh_v1.*kforward.*PZfactor*dt;
    
    ponTh_v2 = ponTh_v2 + PON_v1.*dTh_v1.*kforward*dt;
    dTh_v2 = dTh_v2 - PON_v1.*dTh_v1.*kforward*dt;
    
    lponTh_v2 = lponTh_v2 + LPON_v1.*dTh_v1.*kforward*dt;
    dTh_v2 = dTh_v2 - LPON_v1.*dTh_v1.*kforward*dt;
    
    %------------Decay----------------------------------------------------------------------------------------------
    dTh_v2 = dTh_v2 - dTh_v1.*lam234*dt;
    
    spTh_v2 = spTh_v2 - spTh_v1.*lam234*dt;
    
    lpTh_v2 = lpTh_v2 - lpTh_v1.*lam234*dt;
    
    szTh_v2 = szTh_v2 - szTh_v1.*lam234*dt;
    
    lzresTh_v2 = lzresTh_v2 - lzresTh_v1.*lam234*dt;
    
    lzdvmTh_v2 = lzdvmTh_v2 - lzdvmTh_v1.*lam234*dt;
    
    pzresTh_v2 = pzresTh_v2 - pzresTh_v1.*lam234*dt;
    
    pzdvmTh_v2 = pzdvmTh_v2 - pzdvmTh_v1.*lam234*dt;
    
    ponTh_v2 = ponTh_v2 - ponTh_v1.*lam234*dt;
    
    lponTh_v2 = lponTh_v2 - lponTh_v1.*lam234*dt;
    
    %------------Production----------------------------------------------------------------------------------------------
    U238 = salinity.*0.0786-0.315;            %From Owens et al. (2011)
    dTh_v2 = dTh_v2 + U238.*lam234*dt;
    
end

if N15on==1
    %----N15 Parameters---------------------------------------------
    %Nitrogen Assimilation
    Eps_NO3up=-5;           %
    Eps_NH4up=-10;          %    
    %Phytoplankton Exudation
    Eps_Exu=0;
    %Excretion
    Eps_ZsExc=-1;           %
    Eps_ZlExc=-5;           %
    Eps_ZpExc=-5;          %    
    %Egestion
    Eps_ZsEg=-1;           %
    Eps_ZlEg=-2;           %
    Eps_ZpEg=-2;           %    
    %Nitrification
    Eps_Nit=-14;          %
    %Remineralization
    Eps_P2D=-1;          %
    Eps_P2N=-1;          %
    Eps_D2N=-1;          %
    Eps_D2R=0;          %
    %Nitrogen Isotopic Ratio of N2
    RN2=0.0036765;           %
    
    Alpha_NO3up=exp(Eps_NO3up/1000);
    Alpha_NH4up=exp(Eps_NH4up/1000);
    Alpha_Exu=exp(Eps_Exu/1000);
    Alpha_ZsExc=exp(Eps_ZsExc/1000);
    Alpha_ZlExc=exp(Eps_ZlExc/1000);
    Alpha_ZpExc=exp(Eps_ZpExc/1000);
    Alpha_ZsEg=exp(Eps_ZsEg/1000);
    Alpha_ZlEg=exp(Eps_ZlEg/1000);
    Alpha_ZpEg=exp(Eps_ZpEg/1000);
    Alpha_Nit=exp(Eps_Nit/1000);
    Alpha_P2D=exp(Eps_P2D/1000);
    Alpha_P2N=exp(Eps_P2N/1000);
    Alpha_D2N=exp(Eps_D2N/1000);
    Alpha_D2R=exp(Eps_D2R/1000);
    
    %-----Initialize Variables--------------------------------------
    spN15_v1=tracer(:,iN15+0);       spN15_v2=spN15_v1; 
    lpN15_v1=tracer(:,iN15+1);       lpN15_v2=lpN15_v1;
    szN15_v1=tracer(:,iN15+2);       szN15_v2=szN15_v1;
    lzresN15_v1=tracer(:,iN15+3);    lzresN15_v2=lzresN15_v1;
    lzdvmN15_v1=tracer(:,iN15+4);    lzdvmN15_v2=lzdvmN15_v1;
    pzresN15_v1=tracer(:,iN15+5);    pzresN15_v2=pzresN15_v1;
    pzdvmN15_v1=tracer(:,iN15+6);    pzdvmN15_v2=pzdvmN15_v1;
    NON15_v1=tracer(:,iN15+7);       NON15_v2=NON15_v1;
    NHN15_v1=tracer(:,iN15+8);       NHN15_v2=NHN15_v1;
    PONN15_v1=tracer(:,iN15+9);      PONN15_v2=PONN15_v1;
    LPONN15_v1=tracer(:,iN15+10);     LPONN15_v2=LPONN15_v1;
    DONN15_v1=tracer(:,iN15+11);     DONN15_v2=DONN15_v1;
    DONrefN15_v1=tracer(:,iN15+12);  DONrefN15_v2=DONrefN15_v1;
    
    %-----N15 Ratios in state variables----------------------------
    Rsp = spN15_v1./sp_v1;
    Rlp = lpN15_v1./lp_v1;
    Rsz = szN15_v1./sz_v1;
    Rlzres = lzresN15_v1./lzres_v1;
    Rlzdvm = lzdvmN15_v1./lzdvm_v1;
    Rpzres = pzresN15_v1./pzres_v1;
    Rpzdvm = pzdvmN15_v1./pzdvm_v1;
    RNO3 = NON15_v1./NO_v1;
    RNH4 = NHN15_v1./NH_v1;
    RPON = PONN15_v1./PON_v1;
    RLPON = LPONN15_v1./LPON_v1;
    RDON = DONN15_v1./DON_v1;
    RDONref = DONrefN15_v1./DONref_v1;
    
    %-----Nitrate Uptake-------------------------------------------
    spN15_v2 = spN15_v2 + NO3up_sp.*RNO3.*Alpha_NO3up ...
                        + NH4up_sp.*RNH4.*Alpha_NH4up ...
                        - exu_sp.*Rsp.*Alpha_Exu;
    NON15_v2 = NON15_v2 - NO3up_sp.*RNO3.*Alpha_NO3up;
    NHN15_v2 = NHN15_v2 - NH4up_sp.*RNH4.*Alpha_NH4up;
    DONN15_v2 = DONN15_v2 + exu_sp.*Rsp.*Alpha_Exu;
    
    lpN15_v2 = lpN15_v2 + NO3up_lp.*RNO3.*Alpha_NO3up ...
                        + NH4up_lp.*RNH4.*Alpha_NH4up ...
                        - exu_lp.*Rlp.*Alpha_Exu;
    NON15_v2 = NON15_v2 - NO3up_lp.*RNO3.*Alpha_NO3up;
    NHN15_v2 = NHN15_v2 - NH4up_lp.*RNH4.*Alpha_NH4up;
    DONN15_v2 = DONN15_v2 + exu_lp.*Rlp.*Alpha_Exu;
    
    %-----Phytoplankton respiration and mortality-------------------------------------------
    spN15_v2 = spN15_v2 - respiration_sp.*np_frac_sp.*RNO3.*Alpha_NO3up ...
                        - respiration_sp.*(1-np_frac_sp).*RNH4.*Alpha_NH4up;
    NON15_v2 = NON15_v2 + respiration_sp.*np_frac_sp.*RNO3.*Alpha_NO3up;
    NHN15_v2 = NHN15_v2 + respiration_sp.*(1-np_frac_sp).*RNH4.*Alpha_NH4up;
    
    lpN15_v2 = lpN15_v2 - respiration_lp.*np_frac_lp.*RNO3.*Alpha_NO3up ...
                        - respiration_lp.*(1-np_frac_lp).*RNH4.*Alpha_NH4up;
    NON15_v2 = NON15_v2 + respiration_lp.*np_frac_lp.*RNO3.*Alpha_NO3up;
    NHN15_v2 = NHN15_v2 + respiration_lp.*(1-np_frac_lp).*RNH4.*Alpha_NH4up;
    
    spN15_v2 = spN15_v2 - mortality_sp.*Rsp;
    PONN15_v2 = PONN15_v2 + mortality_sp.*Rsp;
    
    lpN15_v2 = lpN15_v2 - mortality_lp.*Rlp;
    PONN15_v2 = PONN15_v2 + mortality_lp.*Rlp;
    
    %-----Zooplankton Grazing----------------------------------------------------------------
    spN15_v2 = spN15_v2 - SP2SZ.*Rsp;
    lpN15_v2 = lpN15_v2 - LP2SZ.*Rlp;
    szN15_v2 = szN15_v2 + SP2SZ.*Rsp + LP2SZ.*Rlp;
    PONN15_v2 = PONN15_v2 + (1-ae_sz)*(SP2SZ+LP2SZ).*Rsz*Alpha_ZsEg;    %----------------FIX THIS IN OTHER CODES!!  FIX IT!! FIX IT!!!
    szN15_v2 = szN15_v2 - (1-ae_sz)*(SP2SZ+LP2SZ).*Rsz*Alpha_ZsEg;
    NHN15_v2 = NHN15_v2 + (ae_sz-gge_sz)*(SP2SZ+LP2SZ).*Rsz*Alpha_ZsExc;
    szN15_v2 = szN15_v2 - (ae_sz-gge_sz)*(SP2SZ+LP2SZ).*Rsz*Alpha_ZsExc;
    
    spN15_v2 = spN15_v2 - SP2LZres.*Rsp;
    lpN15_v2 = lpN15_v2 - LP2LZres.*Rlp;
    szN15_v2 = szN15_v2 - SZ2LZres.*Rsz;
    lzresN15_v2 = lzresN15_v2 + SP2LZres.*Rsp + LP2LZres.*Rlp + SZ2LZres.*Rsz;
    LPONN15_v2 = LPONN15_v2 + (1-ae_lzres)*(SP2LZres+LP2LZres+SZ2LZres).*Rlzres*Alpha_ZlEg;
    lzresN15_v2 = lzresN15_v2 - (1-ae_lzres)*(SP2LZres+LP2LZres+SZ2LZres).*Rlzres*Alpha_ZlEg;
    NHN15_v2 = NHN15_v2 + (act_res_lzres)*(SP2LZres+LP2LZres+SZ2LZres).*Rlzres*Alpha_ZlExc;
    lzresN15_v2 = lzresN15_v2 - (act_res_lzres)*(SP2LZres+LP2LZres+SZ2LZres).*Rlzres*Alpha_ZlExc;
    
    spN15_v2 = spN15_v2 - SP2LZdvm.*Rsp;
    lpN15_v2 = lpN15_v2 - LP2LZdvm.*Rlp;
    szN15_v2 = szN15_v2 - SZ2LZdvm.*Rsz;
    lzdvmN15_v2 = lzdvmN15_v2 + SP2LZdvm.*Rsp + LP2LZdvm.*Rlp + SZ2LZdvm.*Rsz;
    LPONN15_v2 = LPONN15_v2 + (1-ae_lzdvm)*(SP2LZdvm+LP2LZdvm+SZ2LZdvm).*Rlzdvm*Alpha_ZlEg;
    lzdvmN15_v2 = lzdvmN15_v2 - (1-ae_lzdvm)*(SP2LZdvm+LP2LZdvm+SZ2LZdvm).*Rlzdvm*Alpha_ZlEg;
    NHN15_v2 = NHN15_v2 + (act_res_lzdvm)*(SP2LZdvm+LP2LZdvm+SZ2LZdvm).*Rlzdvm*Alpha_ZlExc;
    lzdvmN15_v2 = lzdvmN15_v2 - (act_res_lzdvm)*(SP2LZdvm+LP2LZdvm+SZ2LZdvm).*Rlzdvm*Alpha_ZlExc;
    
    spN15_v2 = spN15_v2 - SP2PZres.*Rsp;
    lpN15_v2 = lpN15_v2 - LP2PZres.*Rlp;
    szN15_v2 = szN15_v2 - SZ2PZres.*Rsz;
    lzresN15_v2 = lzresN15_v2 - LZres2PZres.*Rlzres;
    lzdvmN15_v2 = lzdvmN15_v2 - LZdvm2PZres.*Rlzdvm;
    pzresN15_v2 = pzresN15_v2 + SP2PZres.*Rsp + LP2PZres.*Rlp + SZ2PZres.*Rsz + ...
                  LZres2PZres.*Rlzres + LZdvm2PZres.*Rlzdvm;
    LPONN15_v2 = LPONN15_v2 + (1-ae_pzres)*(SP2PZres+LP2PZres+SZ2PZres+LZres2PZres+LZdvm2PZres).*Rpzres*Alpha_ZpEg;
    pzresN15_v2 = pzresN15_v2 - (1-ae_pzres)*(SP2PZres+LP2PZres+SZ2PZres+LZres2PZres+LZdvm2PZres).*Rpzres*Alpha_ZpEg;
    NHN15_v2 = NHN15_v2 + (act_res_pzres)*(SP2PZres+LP2PZres+SZ2PZres+LZres2PZres+LZdvm2PZres).*Rpzres*Alpha_ZpExc;
    pzresN15_v2 = pzresN15_v2 - (act_res_pzres)*(SP2PZres+LP2PZres+SZ2PZres+LZres2PZres+LZdvm2PZres).*Rpzres*Alpha_ZpExc;
    
    spN15_v2 = spN15_v2 - SP2PZdvm.*Rsp;
    lpN15_v2 = lpN15_v2 - LP2PZdvm.*Rlp;
    szN15_v2 = szN15_v2 - SZ2PZdvm.*Rsz;
    lzresN15_v2 = lzresN15_v2 - LZres2PZdvm.*Rlzres;
    lzdvmN15_v2 = lzdvmN15_v2 - LZdvm2PZdvm.*Rlzdvm;
    pzdvmN15_v2 = pzdvmN15_v2 + SP2PZdvm.*Rsp + LP2PZdvm.*Rlp + SZ2PZdvm.*Rsz + ...
                  LZres2PZdvm.*Rlzres + LZdvm2PZdvm.*Rlzdvm;
    LPONN15_v2 = LPONN15_v2 + (1-ae_pzdvm)*(SP2PZdvm+LP2PZdvm+SZ2PZdvm+LZres2PZdvm+LZdvm2PZdvm).*Rpzdvm*Alpha_ZpEg;
    pzdvmN15_v2 = pzdvmN15_v2 - (1-ae_pzdvm)*(SP2PZdvm+LP2PZdvm+SZ2PZdvm+LZres2PZdvm+LZdvm2PZdvm).*Rpzdvm*Alpha_ZpEg;
    NHN15_v2 = NHN15_v2 + (act_res_pzdvm)*(SP2PZdvm+LP2PZdvm+SZ2PZdvm+LZres2PZdvm+LZdvm2PZdvm).*Rpzdvm*Alpha_ZpExc;
    pzdvmN15_v2 = pzdvmN15_v2 - (act_res_pzdvm)*(SP2PZdvm+LP2PZdvm+SZ2PZdvm+LZres2PZdvm+LZdvm2PZdvm).*Rpzdvm*Alpha_ZpExc;
    
    %-----Zooplankton Mortality----------------------------------------------------------------
    szN15_v2 = szN15_v2 - mortality_sz.*Rsz;
    lzresN15_v2 = lzresN15_v2 - mortality_lzres.*Rlzres;
    lzdvmN15_v2 = lzdvmN15_v2 - mortality_lzdvm.*Rlzdvm;
    pzresN15_v2 = pzresN15_v2 - mortality_pzres.*Rpzres;
    pzdvmN15_v2 = pzdvmN15_v2 - mortality_pzdvm.*Rpzdvm;
    if day==1 & euponly==1
        LPONN15_v2 = LPONN15_v2 + mortality_sz.*Rsz + mortality_lzres.*Rlzres + mortality_pzres.*Rpzres;
    else
        LPONN15_v2 = LPONN15_v2 + mortality_sz.*Rsz + mortality_lzres.*Rlzres + mortality_lzdvm.*Rlzdvm ...
                                + mortality_pzres.*Rpzres + mortality_pzdvm.*Rpzdvm;
    end
    
    %-----Zooplankton Basal Excretion----------------------------------------------------------------
    lzresN15_v2 = lzresN15_v2 - basalexcretion_lzres.*Rlzres*Alpha_ZlExc;
    lzdvmN15_v2 = lzdvmN15_v2 - basalexcretion_lzdvm.*Rlzdvm*Alpha_ZlExc;
    pzresN15_v2 = pzresN15_v2 - basalexcretion_pzres.*Rpzres*Alpha_ZpExc;
    pzdvmN15_v2 = pzdvmN15_v2 - basalexcretion_pzdvm.*Rpzdvm*Alpha_ZpExc;
    if day==1 & euponly==1
        NHN15_v2 = NHN15_v2 + basalexcretion_lzres.*Rlzres*Alpha_ZlExc + basalexcretion_pzres.*Rpzres*Alpha_ZpExc;
    else
        NHN15_v2 = NHN15_v2 + basalexcretion_lzres.*Rlzres*Alpha_ZlExc + basalexcretion_lzdvm.*Rlzdvm*Alpha_ZlExc ...
                            + basalexcretion_pzres.*Rpzres*Alpha_ZpExc + basalexcretion_pzdvm.*Rpzdvm*Alpha_ZpExc;
    end
    
    %-----Chemical Transformations----------------------------------------------------------------
    NHN15_v2 = NHN15_v2 - nitrification_NH.*RNH4*Alpha_Nit;
    NON15_v2 = NON15_v2 + nitrification_NH.*RNH4*Alpha_Nit;
    
    PONN15_v2 = PONN15_v2 - aggregation_PON.*RPON;
    LPONN15_v2 = LPONN15_v2 + aggregation_PON.*RPON;
    
    PONN15_v2 = PONN15_v2-remineralization_PON.*RPON*Alpha_P2N;
    NHN15_v2 = NHN15_v2+remineralization_PON.*RPON*Alpha_P2N;
    
    LPONN15_v2 = LPONN15_v2-remineralization_LPON.*RLPON*Alpha_P2N;
    NHN15_v2 = NHN15_v2+remineralization_LPON.*RLPON*Alpha_P2N;
    
    DONN15_v2 = DONN15_v2-remineralization_DON*(1-f_DON_DONref).*RDON*Alpha_D2N;
    NHN15_v2 = NHN15_v2+remineralization_DON*(1-f_DON_DONref).*RDON*Alpha_D2N;
    DONN15_v2 = DONN15_v2-remineralization_DON*f_DON_DONref.*RDON;   %No fractionation of the portion of DON that is left behind as refractory DON
    DONrefN15_v2 = DONrefN15_v1+remineralization_DON*f_DON_DONref.*RDON;     %No fractionation of the portion of DON that is left behind as refractory DON
    
%     DONN15_v2 = DONN15_v2-DON2Refractory.*RDON*Alpha_D2R;
%     DONrefN15_v2 = DONrefN15_v1+DON2Refractory.*RDON*Alpha_D2R;

    DONrefN15_v2 = DONrefN15_v2-remineralization_DONref.*RDONref*Alpha_D2N;
    NHN15_v2 = NHN15_v2+remineralization_DONref.*RDONref*Alpha_D2N;
    
    PONN15_v2 = PONN15_v2-dissolution_PON.*RPON*Alpha_P2D;
    DONN15_v2 = DONN15_v2+dissolution_PON.*RPON*Alpha_P2D;
    
    LPONN15_v2 = LPONN15_v2-dissolution_LPON.*RLPON*Alpha_P2D;
    DONN15_v2 = DONN15_v2+dissolution_LPON.*RLPON*Alpha_P2D;
    
end
    
    

tracer_out(:,1)=sp_v2;
tracer_out(:,2)=lp_v2;
tracer_out(:,3)=sz_v2;
tracer_out(:,4)=lzres_v2;
tracer_out(:,5)=lzdvm_v2;
tracer_out(:,6)=pzres_v2;
tracer_out(:,7)=pzdvm_v2;
tracer_out(:,8)=NO_v2;
tracer_out(:,9)=NH_v2;
tracer_out(:,10)=PON_v2;
tracer_out(:,11)=LPON_v2;
tracer_out(:,12)=DON_v2;
tracer_out(:,13)=DONref_v2;
tracer_out(:,14)=SI_v2;
tracer_out(:,15)=OP_v2;
tracer_out(:,16)=LOP_v2;
tracer_out(:,17)=Chl_ps_v2;
tracer_out(:,18)=Chl_pl_v2;
tracer_out(:,19)=oxy_v2;
if CO2on==1
    tracer_out(:,iCO2+0)=DIC_v2;
    tracer_out(:,iCO2+1)=alk_v2;
end
if Thoriumon==1
    tracer_out(:,iThorium+0)=dTh_v2;
    tracer_out(:,iThorium+1)=spTh_v2;
    tracer_out(:,iThorium+2)=lpTh_v2;
    tracer_out(:,iThorium+3)=szTh_v2;
    tracer_out(:,iThorium+4)=lzresTh_v2;
    tracer_out(:,iThorium+5)=lzdvmTh_v2;
    tracer_out(:,iThorium+6)=pzresTh_v2;
    tracer_out(:,iThorium+7)=pzdvmTh_v2;
    tracer_out(:,iThorium+8)=ponTh_v2;
    tracer_out(:,iThorium+9)=lponTh_v2;
end
if N15on==1
    tracer_out(:,iN15+0)=       spN15_v2; 
    tracer_out(:,iN15+1)=       lpN15_v2;
    tracer_out(:,iN15+2)=       szN15_v2;
    tracer_out(:,iN15+3)=    lzresN15_v2;
    tracer_out(:,iN15+4)=    lzdvmN15_v2;
    tracer_out(:,iN15+5)=    pzresN15_v2;
    tracer_out(:,iN15+6)=    pzdvmN15_v2;
    tracer_out(:,iN15+7)=       NON15_v2;
    tracer_out(:,iN15+8)=       NHN15_v2;
    tracer_out(:,iN15+9)=      PONN15_v2;
    tracer_out(:,iN15+10)=     LPONN15_v2;
    tracer_out(:,iN15+11)=     DONN15_v2;
    tracer_out(:,iN15+12)=  DONrefN15_v2;
end

basalexcretion_dvm = basalexcretion_lzdvm + basalexcretion_pzdvm;
activexcretion_dvm = (act_res_lzdvm)*(SP2LZdvm+LP2LZdvm+SZ2LZdvm) + (act_res_pzdvm)*(SP2PZdvm+LP2PZdvm+SZ2PZdvm+LZres2PZdvm+LZdvm2PZdvm);
mortality_dvm = mortality_lzdvm + mortality_pzdvm;