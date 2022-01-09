function [Param]=InitParameters(input)

% Small Phytoplankton
Param(1)=0.0693;  %tdep = Temperature Dependence
Param(2)=0.4;     %vmax_sp = max growth rate at 0 deg C
Param(3)=0.5;     %k_NO_sp = 0.5;                  % nitrate half saturation constant - [1.0]
Param(4)=0.1;     %k_NH_sp = 0.1;                  % ammonium half saturation constant - [0.1]
Param(5)=0.1;     %alpha_sp = 0.1;    	    	% PI curve parameter alpha - [0.01]
Param(6)=4.5e-4;   %beta_sp = 4.5e-4;                 % PI curve parmater beta - [4.5e-4], 1.4e-3
Param(7)=0.03;     %ref_resp_sp = 0.03;            % respiration at 0 deg C (1/d) - [0.03]
Param(8)=0.002;   %ref_mort_sp = 0.002;              % mortality at at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
Param(9)=0.135;   %ext_excr_sp = 0.135;             % extracellular excretion - [0.135]
Param(10)=1.4;    %inh_NH_NO_sp = 1.4;             %Ammonium inhibition of nitrate uptake

% Large Phytoplankton 
Param(11)=0.8;    %vmax_lp = 0.8;                  % max growth rate at 0 deg C (1/d) - [0.8]
Param(12)=3;      %k_NO_lp = 3.0;                  % nitrate half saturation constant - [3.0]
Param(13)=0.3;    %k_NH_lp = 0.3;                  % ammonium half saturation constant - [0.3]
Param(14)=6;      %k_SI_lp = 6.0;                  % silica half saturation constant - [6.0] 
Param(15)=0.006;    %alpha_lp = 0.1;                 % PI curve parameter alpha - [0.01], 1.4e-3
Param(16)=4.5e-4; %beta_lp= 4.5e-4;                  % PI curve parmater beta - [4.5e-4]
Param(17)=0.03;   %ref_resp_lp = 0.03;             % respiration at 0 deg C (1/d) - [0.03]
Param(18)=0.001;  %ref_mort_lp = 0.001;              % mortality at at 0 deg C (1/d) - [0.029], Linear Mort - ~[0.001]
Param(19)=0.135;  %ext_excr_lp = 0.135;              % extracellular excretion - [0.0135]
Param(20)=1.4;    %inh_NH_NO_lp = 1.4;              %Ammonium inhibition of nitrate utpake

% Small Zooplankton 
Param(21)=0.6;     %gmax_sz_sp = 0.6;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.04]
Param(22)=0.6;     %gmax_sz_lp = 0.6;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.04]
Param(23)=1.4;     %iv_sz_sp = 1.4;                 % ivlev constant - [1.4]
Param(24)=1.4;     %iv_sz_lp = 1.4;                 % ivlev constant - [1.4]
Param(25)=0.04;    %thresh_sz_sp = 0.04 ;           % feeding threshold on small phytoplankton - [0.043]
Param(26)=0.04;    %thresh_sz_lp = 0.04 ;           % feeding threshold on large phytoplankton - [0.043]
Param(27)=0.022;   %ref_mort_sz = 0.022;              % mortality at 0 deg C (1/d) - [0.0585], Linear Mort - ~[0.002]
Param(28)=0.70;    %ae_sz = 0.70;                   % assimilation efficency - [0.7]
Param(29)=0.3;     %gge_sz = 0.30;                  % gross growth efficency - [0.3]

% Large Zooplankton - Resident
Param(30)=0.1;       %gmax_lzres_sp = 0.1;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.1]
Param(31)=0.3;     %gmax_lzres_lp = 0.3;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.4]
Param(32)=0.3;     %gmax_lzres_sz = 0.3;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.4]
Param(33)=1.4;     %iv_lz_sp = 1.4;                 % ivlev constant - [1.4]
Param(34)=1.4;     %iv_lz_lp = 1.4;                 % ivlev constant - [1.4]
Param(35)=1.4;     %iv_lz_sz = 1.4;                 % ivlev constant - [1.4]
Param(36)=0.7;     %ae_lzres = 0.30;                  % assimilation efficiency of lzres
Param(37)=0.1;     %act_res_lzres = 0.1;            %Active respiration of lzres
Param(38)=0.1;     %mort_day_lzres                  %Daytime mortality parameter for lzres

Param(39)=0.04;    %thresh_lz = 0.04;            % feeding threshold on small phytoplankton - [0.04]
Param(40)=0.04;    %thresh_lz = 0.04;            % feeding threshold on large phytoplankton - [0.04]
Param(41)=0.04;    %thresh_lz = 0.04;            % feeding threshold on small zooplankton  - [0.04]
Param(42)=0.1;     %mort_night_lz                  %Nighttime mortality parameter for lzres
Param(43)=0.0656;     %Ikeda_a2
Param(44)=0.0296;     %Ikeda_lz

% Large Zooplankton - Migratory
Param(45)=0.1;       %gmax_lzdvm_sp = 0.1;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.1]
Param(46)=0.3;     %gmax_lzdvm_lp = 0.3;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.4]
Param(47)=0.3;     %gmax_lzdvm_sz = 0.3;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.4]
Param(48)=1.4;     %iv_lz_sp = 1.4;                 % ivlev constant - [1.4]
Param(49)=1.4;     %iv_lz_lp = 1.4;                 % ivlev constant - [1.4]
Param(50)=1.4;     %iv_lz_sz = 1.4;                 % ivlev constant - [1.4]
Param(51)=0.7;     %ae_lzres = 0.30;                  % assimilation efficiency of lzres
Param(52)=0.1;     %act_res_lzres = 0.1;            %Active respiration of lzres
Param(53)=0.1;     %mort_day_lzres                  %Daytime mortality parameter for lzres

% Predatory Zooplankton - Resident
Param(54)=0.1;     %gmax_pz_sp = 0.1;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.2]
Param(55)=0.1;     %gmax_pz_lp = 0.1;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.2]
Param(56)=0.1;     %gmax_pz_sz = 0.1;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.2]
Param(57)=0.3;     %gmax_pz_lz = 0.3;               % max grazing rate on large zooplankton at 0 deg C (1/d) - [0.2]
Param(58)=1.4;     %iv_pz_sp = 1.4;                 % ivlev constant - [1.4]
Param(59)=1.4;     %iv_pz_lp = 1.4;                 % ivlev constant - [1.4]
Param(60)=1.4;     %iv_pz_sz = 1.4;                 % ivlev constant - [1.4]
Param(61)=1.4;     %iv_pz_lz = 1.4;                 % ivlev constant - [1.4]
Param(62)=0.7;     %ae_pz = 0.70;                    % assimilation efficency - [0.7]
Param(63)=0.1;     %act_res_pzres;                  % active respiration (percent of grazing)
Param(64)=0.1;     %mort_day_pzres                  %Nighttime mortality parameter for lzres

Param(65)=4.605;   %inh_szlzlp_sp = 4.605;            % grazing on large phytoplankton inhibition by small and large zooplankton - [4.605]
Param(66)=4.605;   %inh_szlz_lp = 4.605;            % grazing on large phytoplankton inhibition by small and large zooplankton - [4.605]
Param(67)=3.01;    %inh_lz_sz = 3.01;               % grazing on small zooplankton inhibition by large zooplankton - [3.01]
Param(68)=0.04;    %thresh_pz_lp = 0.04;            % feeding threshold on small phytoplankton - [0.04]
Param(69)=0.04;    %thresh_pz_lp = 0.04;            % feeding threshold on large phytoplankton - [0.04]
Param(70)=0.04;    %thresh_pz_sz = 0.04;            % feeding threshold on small zooplankton - [0.04]
Param(71)=0.04;    %thresh_pz_lz = 0.04;            % feeding threshold on large zooplankton - [0.04]
Param(72)=0.05;     %mort_night_pz                  %Nighttime mortality parameter for lzres
Param(73)=0.0173;                   %Ikeda_pz

% Predatory Zooplankton - Migratory
Param(74)=0.2;     %gmax_pzdvm_sp = 0.1;               % max grazing rate on small phytoplankton at 0 deg C (1/d) - [0.2]
Param(75)=0.2;     %gmax_pzdvm_lp = 0.1;               % max grazing rate on large phytoplankton at 0 deg C (1/d) - [0.2]
Param(76)=0.2;     %gmax_pzdvm_sz = 0.1;               % max grazing rate on small zooplankton at 0 deg C (1/d) - [0.2]
Param(77)=0.3;     %gmax_pzdvm_lz = 0.3;               % max grazing rate on large zooplankton at 0 deg C (1/d) - [0.2]
Param(78)=1.4;     %iv_pzdvm_sp = 1.4;                 % ivlev constant - [1.4]
Param(79)=1.4;     %iv_pzdvm_lp = 1.4;                 % ivlev constant - [1.4]
Param(80)=1.4;     %iv_pzdvm_sz = 1.4;                 % ivlev constant - [1.4]
Param(81)=1.4;     %iv_pzdvm_lz = 1.4;                 % ivlev constant - [1.4]
Param(82)=0.7;     %ae_pzdvm = 0.70;                    % assimilation efficency - [0.7]
Param(83)=0.1;     %act_res_pzdvm;                  % active respiration (percent of grazing)
Param(84)=0.05;     %mort_day_pzdvm                  %Nighttime mortality parameter for lzres

% Nutrients (11)
Param(85)=0.003;   %ref_nitr = 0.003;                % nitrification at 0 deg C (1/d) - [0.03]
Param(86)=0.01;    %ref_dec_PON_NH = 0.01;           % decompositon from PON to NH at 0 deg C (1/d) - [0.1]
Param(87)=0.05;    %ref_dec_PON_DON = 0.05;          % decompositon from PON to DON at 0 deg C (1/d) - [0.1]
Param(88)=0.05;    %ref_dec_PON_LPON = 0.05;          % aggregation from PON to LPON at 0 deg C (1/d) - [0.1]
Param(89)=0.01;    %ref_dec_LPON_NH = 0.01;           % decompositon from LPON to NH at 0 deg C (1/d) - [0.1]
Param(90)=0.05;    %ref_dec_LPON_DON = 0.05;          % decompositon from LPON to DON at 0 deg C (1/d) - [0.1]
Param(91)=0.01;    %ref_dec_DON_NH = 0.02;           % decompositon from DON to NH at 0 deg C (1/d) - [0.2] 
Param(92)=0.01;    %ref_dec_DON_DONref = 0.02;           % decompositon from DON to NH at 0 deg C (1/d) - [0.2] 
Param(93)=10^-6;    %ref_dec_DONref_NH = 0.02;           % decompositon from DON to NH at 0 deg C (1/d) - [0.2] 
Param(94)=0.01;    %ref_dec_OP_SI = 0.01;            % decompositon from OP to SI at 0 deg C (1/d) - [0.1]
Param(95)=15;      %sink=15;                    		% sinking speed per day - [40] 
Param(96)=100;      %Lsink=15;                    		% sinking speed per day of large PON and OP - [40] 

% Oxygen Limitation
Param(97) = 20;     %koxy_mic = 20 umol / kg           %Oxygen limitation (microbes = protists + implicit)
Param(98) = 100;    %koxy_met = 100 umol / kg          %Oxygen limitation (metazoans)

% Chl2c Submodel (Li et al., 2010)
Param(99)=0.001;         %chl2c_sp_min = 0.0001; 		% minimum Chl:C ratio - [0.000]
Param(100)=0.005;         %chl2c_lp_min = 0.005; 		% maximum Chl:C ratio - [0.005] 
Param(101)=0.03;         %chl2c_sp_max = 0.03; 		% minimum Chl:C ratio - [0.03]
Param(102)=0.065;         %chl2c_lp_max = 0.061; 		% maximum Chl:C ratio - [0.061] 

%Light
Param(103)=0.05;         %ext_chl                          %Light absorption by Chl (m-1 / (mg Chl m-3)) - 0.05



