function [ddo] = oxy_sat(oxy_v1,temperature,salinity)


A0=2.00907;
A1= 3.22014;
A2= 4.05010;
A3= 4.94457;
A4=-2.56847*10^-1;
A5= 3.88767 ;

B0=-6.24523*10^-3;
bB1=-7.37614*10^-3;
B2=-1.03410*10^-2;
B3=-8.17083*10^-3;
C0=-4.88682*10^-7;

TT  = 298.15-temperature;
TK  = 273.15+temperature;
S = salinity;
TS  = log(TT/TK);
TS2 = TS^2;
TS3 = TS^3;
TS4 = TS^4;
TS5 = TS^5;
CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5 ...
    + S*(B0 + bB1*TS + B2*TS2 + B3*TS3) ...
    + C0*(S*S);
o2sato = exp(CO);
o2sato = (o2sato*43.6*1025/(1000*1000));
ddo = o2sato - (oxy_v1/(1000.0));