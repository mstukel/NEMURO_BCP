function [tracer_out,co2flux,surfpco2,diffpco2] = Carbon(tracer,temperature,salinity,atmco2,thickness,ph,wind,dt)

%Note that dt should be in units of seconds
dt = dt*24*60*60;

atmpres = 1;

DIC_v1=tracer(1,20)/1000;           %Dissolved inorganic carbon (units = mol / m3) in surface layer
alk_v1=tracer(1,21)/1000;           %Alkalinity (units = mol / m3) in surface layer
SI_v1 =tracer(1,14)/1000;           %Silicic acid (units = mol / m3) in surface layer
oxy_v1 = tracer(1,19);              %Oygen
PO4_est = (tracer(1,8)+tracer(1,9))/16/1000;   %estimate phosphate from NO3 + NH4


if ph<7 | ph>9.5
    ph = 8;
end
atmospco2 = atmco2*10^-6;


pressure=0;



[co2star,dco2star,pco2surf,dpco2] = co2cal_SWS(DIC_v1,alk_v1,SI_v1,PO4_est,pressure,temperature,salinity,ph,atmospco2,atmpres);
diff_co2 = dco2star;
[diff_o2] = oxy_sat(oxy_v1,temperature,salinity);
[gflux_co2,gflux_oxy,tv_oxy,tv_co2] = gasexch(wind,temperature,thickness,diff_co2,diff_o2); %Units are moles/m2/second

DIC_v2 = DIC_v1 + gflux_co2*dt;
oxy_v2 = oxy_v1 + gflux_oxy*1000*dt;
surfpco2 = pco2surf * 10^6;
diffpco2 = dpco2 * 10^6;
co2flux = gflux_co2 * 60*60*24*365*12;  %grams / m2 / yr

tracer_out = tracer;
tracer_out(1,20)=DIC_v2*1000;
%tracer_out(1,21)=alk_v2*1000;
tracer_out(1,19)=oxy_v2;


