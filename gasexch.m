function [gflux_co2,gflux_oxy,tv_oxy,tv_co2] = gasexch(wind,temperature,thickness,diff_co2,diff_o2)


scc = 660;
sco = 663;
bc = 49;
bo = 9;

%-------CO2
sc = 2073.1-125.62*temperature+ ...
    (3.6276*temperature*temperature) - ...
    (0.043219*temperature*temperature*temperature);
sc = (sc/scc);
skin = (-0.0005*temperature + 0.008);
bubble = 0.01*(wind/bc)*(wind/bc);

%Schmidt number for piston or gas transfer velocity (cm/hr)
%cm/hr to m/sec = 1/(100*60*60) = .0000027777
tv_co2 = 0.31*wind*wind*(sc^(-1/2))*0.0000027777;
gflux_co2 = tv_co2*(1+skin+bubble)*(diff_co2)/thickness;



%-------oxygen

sc = 1638.-81.83*temperature + ...
    (1.483*temperature*temperature) - ...
    (0.008004*temperature*temperature*temperature);
sc = (sc/sco);
skin = (-0.0004*temperature + 0.0053);
bubble = 0.01*(wind/bo)*(wind/bo);


%Schmidt number for piston or gas transfer velocity (cm/hr)
%cm/hr to m/sec = 1/(100*60*60) = .0000027777
tv_oxy = 0.31*wind*wind*(sc^(-1/2))*0.0000027777;
gflux_oxy = tv_oxy*(1+skin+bubble)*(diff_o2)/thickness;
