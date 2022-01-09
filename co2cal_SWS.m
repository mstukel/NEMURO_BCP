function [co2star,dco2star,pCO2surf,dpCO2,ph] = co2cal_SWS(dic_in,alk_in,SI_in,PO4_in,pressure,temperature,salinity,ph,xco2_in,atmpres)


% --- pressure dependence following MIT code
%-------------------------------------------------------------------------
% NOW; "All" constants are given on seawater H scale (hSWS)
% - HOWEVER, dissociation constants for S and F are on the 'free' H scale
%            (this is necessary to work with hSWS)
% SUBROUTINE CO2CALC_SWS
% PURPOSE
%	Calculate delta co2* from total alkalinity and total CO2 at
% temperature (t), salinity (s) and "atmpres" atmosphere total pressure.

% INPUT
%	dic_in = total inorganic carbon (mol/m^3)
%                where 1 T = 1 metric ton = 1000 kg
%	ta_in  = total alkalinity (eq/m^3)
%	pt_in  = inorganic phosphate (mol/m^3)
%	sit_in = inorganic silicate (mol/m^3)
%	t      = temperature (degrees C)
%	s      = salinity (PSU)
%	phlo   = lower limit of pH range
%	phhi   = upper limit of pH range
%	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) *1.e-6
%	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
% OUTPUT
%	co2star  = CO2*water (mol/m^3)
%	dco2star = delta CO2 (mol/m^3)
%       pco2surf = oceanic pCO2 (ppmv)*1.e-6
%       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)*1.e-6
%--------------------------------------------------------------------------
%

%-- determine pressure from depth use upper surface of cell (1atm = 101325pa)
pressc = pressure/101325.0+1.0;
vomega = 5.0;

%       ---------------------------------------------------------------------
%       Change units from the input of mmol/m^3 -> mol/kg:
%       (1 mmol/m^3)  x (1 m^3/1024.5 kg) x (1mol/1000mmol)
%       where the ocean's mean surface density is 1024.5 kg/m^3
%       Note: mol/kg are actually what the body of this routine uses
%       for calculations.
%       ---------------------------------------------------------------------
%        permil = 1.0 / (1024.5*1.e3)
permil = .000976085;
pt=PO4_in*permil;
sit=SI_in*permil;
ta=alk_in*permil;
dic=dic_in*permil;
%       ---------------------------------------------------------------------
%       Change units from uatm to atm. That is, atm is what the body of
%       this routine uses for calculations.
%       ---------------------------------------------------------------------
permeg=1.e-6;
%       To convert input in uatm -> atm
xco2=xco2_in;
%*************************************************************************
% Calculate all constants needed to convert between various measured
% carbon species. References for each equation are noted in the code.
% Once calculated, the constants are
% stored and passed in the common block "const". The original version of this
% code was based on the code by Dickson in Version 2 of "Handbook of Methods
% for the Analysis of the Various Parameters of the Carbon Dioxide System
% in Seawater", DOE, 1994 (SOP No. 3, p25-26).
%
% Derive simple terms used more than once
%
tk = 273.15 + temperature;
tk100 = tk/100.0;
tk1002=tk100*tk100;
invtk=1.0/tk;
dlogtk=log(tk);
is=19.924*salinity/(1000.-1.005*salinity);
is2=is*is;
sqrtis=sqrt(is);
s2=salinity*salinity;
sqrts=sqrt(salinity);
s15=salinity^1.5;
scl=salinity/1.80655;
%
%------------------------------------------------------------------------
% Calculate concentrations for borate, sulfate, and fluoride
%
% Uppstrom (1974)
bt = 0.000232 * scl/10.811;
% Morris & Riley (1966)
st = 0.14 * scl/96.062;
% Riley (1965)
ft = 0.000067 * scl/18.9984;
%
%------------------------------------------------------------------------
% f = k0(1-pH2O)*correction term for non-ideality
%
% Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
%
ff = exp(-162.8301 + 218.2968/tk100  + ...
    90.9241*log(tk100) - 1.47696*tk1002 + ...
    salinity * (.025695 - .025225*tk100 + ...
    0.0049867*tk1002));
%
% K0 from Weiss 1974
%
k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) + ...
    salinity * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002));

%
%------------------------------------------------------------------------
% k1 = [H][HCO3]/[H2CO3]
% k2 = [H][CO3]/[HCO3]     on hSWS
%
% Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale
% (Original reference: Dickson and Millero, DSR, 1987)
%
k1=10^(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk - ...
    0.0118 * salinity + 0.000116*s2));

k2=10^(-1*(1394.7*invtk + 4.777 - ...
    0.0184*salinity + 0.000118*s2));

if pressure > 0
    % -- pressure dependence
    % -- Takahashi (1981 geosecs report)
    k1=k1*exp((24.2-0.085*temperature)*(pressc-1.0)/(83.143*tk));
    % -- Co2sys docc Lewis and Wallace (1998)
    k2 = k2*exp((16.4-0.040*temperature)*(pressc-1.0)/(83.143*tk));
end

%------------------------------------------------------------------------
% k1p = [H][H2PO4]/[H3PO4] on hSWS
%
% Millero p.670 (1995)
%
k1p = exp(-4576.752*invtk + 115.540 - 18.453 * dlogtk + ...
    (-106.736*invtk + 0.69171) * sqrts + ...
    (-0.65643*invtk - 0.01844) * salinity);
%
%------------------------------------------------------------------------
% k2p = [H][HPO4]/[H2PO4] on hSWS
%
% Millero p.670 (1995)
%
k2p = exp(-8814.715*invtk + 172.1033 - 27.927 * dlogtk + ...
    (-160.340*invtk + 1.3566) * sqrts + ...
    (0.37335*invtk - 0.05778) * salinity);
%
%------------------------------------------------------------------------
% k3p = [H][PO4]/[HPO4] on hSWS
%
% Millero p.670 (1995)
%
k3p = exp(-3070.75*invtk - 18.126 + ...
    (17.27039*invtk + 2.81197) * ...
    sqrts + (-44.99486*invtk - 0.09984) * salinity);

%------------------------------------------------------------------------
% ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
%
% Millero p.671 (1995) using data from Yao and Millero (1995)
% change to (mol/ kg soln)

ksi = exp(-8904.2*invtk + 117.400 - 19.334 * dlogtk + ...
    (-458.79*invtk + 3.5913) * sqrtis + ...
    (188.74*invtk - 1.5998) * is + ...
    (-12.1652*invtk + 0.07871) * is2 + ...
    log(1.0-0.001005*salinity));

%------------------------------------------------------------------------
% kw = [H][OH] on hSWS
%
% Millero p.670 (1995) using composite data
%
kw = exp(-13847.26*invtk + 148.9802 - 23.6521 * dlogtk + ...
    (118.67*invtk - 5.977 + 1.0495 * dlogtk) * ...
    sqrts - 0.01615 * salinity);

%------------------------------------------------------------------------
% ks = [H][SO4]/[HSO4] on free H scale
%
% Dickson (1990, J. chem. Thermodynamics 22, 113)
% change to (mol/ kg soln)

ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk + ...
    (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis + ...
    (35474*invtk - 771.54 + 114.723*dlogtk) * is - ...
    2698*invtk*is^0.5 + 1776*invtk*is2 + ...
    log(1.0 - 0.001005*salinity));

%------------------------------------------------------------------------
% kf = [H][F]/[HF] on free H scale
%
% Dickson and Riley (1979)
% change to (mol/ kg soln)

kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis + ...
    log(1.0 - 0.001005*salinity));

%------------------------------------------------------------------------
% kb = [H][BO2]/[HBO2] on hSWS
%
% Dickson p.673 (1990)
% change from htotal to hSWS

kb=exp( (-8966.90 - 2890.53*sqrts - 77.942*salinity + ...
    1.728*s15 - 0.0996*s2)*invtk + ...
    (148.0248 + 137.1942*sqrts + 1.62142*salinity) + ...
    (-24.4344 - 25.085*sqrts - 0.2474*salinity) * ...
    dlogtk + 0.053105*sqrts*tk + ...
    log((1+(st/ks)+(ft/kf))/(1+(st/ks))) );

% - pressure dependence
% - millero 1995 p. 675
if pressure > 0
    dv = -29.48+0.1622*t+.002608*temperature*temperature;
    pfactor = -(dv/(83.145*tk))*pressc ...
        +(-0.00142/(83.145*tk))*pressc*pressc;
    kb=kb*exp(pfactor);
    
    % -- solubility product for calcite
    % --- mucci(1983)
    
    tmpa1 = -171.9065-(0.077993*tk)+(2839.319/tk) ...
        +(71.595*log10(tk));
    tmpa2 = +(-.77712+(0.0028426*tk)+(178.34/tk))*sqrts;
    tmpa3 = -(0.07711*salinity)+(0.0041249*s15);
    lksp_t_calc=(tmpa1+tmpa2+tmpa3);
    ksp_t_calc=10.0^(lksp_t_calc);
    pressc=(pressc*10.0-10.0)/10.0;
    xvalue=((48.8-0.53*temperature)*pressc ...
        +(-0.00588+0.0001845*temperature)*pressc*pressc) ...
        /(188.93*tk);
    ksp_t_calc=ksp_t_calc*10.0^(xvalue);
end

% --
hSWS = 10.0^(-ph);

%-----------------------------------------------------------------------------------------------------------------------------
[co2star,hSWS] = calc_pco2(pt,sit,ta,dic,hSWS,ff,bt,k1,k2,k1p,k2p,k3p,kb,kw,ksi);
co2starair=ff*xco2*atmpres;
dco2star=co2starair-co2star;


ph=-log10(hSWS);

% -- bicarbonate
if pressure>0
    co3=k1*k2*dic/(hSWS*hSWS+k1*hSWS+k1*k2);
    calcium = .01028*s/35.0;
    vomega = calcium*co3/ksp_t_calc;
end

%       ---------------------------------------------------------------
%      Add two output arguments for storing pCO2surf
%     Should we be using K0 or ff for the solubility here?
%       ---------------------------------------------------------------
pCO2surf = co2star/ff ;
dpCO2    = ((xco2*atmpres)-(pCO2surf));
pCO2surf = co2star/k0 ;
%  Convert units of output arguments
%      Note: co2star and dco2star are calculated in mol/kg within this routine
%      Thus Convert now from mol/kg -> mol/m^3
co2star  = co2star / permil;
dco2star = dco2star / permil;
dic = dic/permil;
ta = ta/permil;
%      Note: pCO2surf and dpCO2 are calculated in atm above.
%      Thus convert now to uatm
%       pCO2surf = pCO2surf / permeg
%       dpCO2    = dpCO2 / permeg