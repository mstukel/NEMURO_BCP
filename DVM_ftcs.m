function[tracernew]=DVM_ftcs(tracer,Coeff0,Coeff1,Coeff2)

tracer(end+1,:)=0;

tracernew = tracer(1:end-1,:).*Coeff0 + ...
    [zeros(1,length(tracer(1,:)));tracer(1:end-2,:)].*Coeff1 + ...
    tracer(2:end,:).*Coeff2;

%profile_new = profile.*Coeff1 + [zeros(1,length(profile(1,:)));profile(1:end-1,:)].*Coeff2

%bottomflux=(tracer(end,:)-tracer(end-1,:))*BottomCoeff;
