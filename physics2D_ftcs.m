function[tracernew,bottomflux]=physics2D_ftcs(tracer,deep,Coeff0,Coeff1,Coeff2,Coeff3,Coeff4)

tracer(end+1,:,:)=deep;


% tracernew = tracer(1:end-1,:,:).*Coeff0 + ...
%     [zeros(size(tracer(1,:,:)));tracer(1:end-2,:,:)].*Coeff1 + ...
%     tracer(2:end,:,:).*Coeff2 + ...
%     [tracer(1:end-1,:,:)] .* Coeff3  + ...
%     [tracer(1:end-1,:,:)] .* Coeff4;

tracernew = tracer(1:end-1,:,:).*Coeff0 + ...
    [zeros(size(tracer(1,:,:)));tracer(1:end-2,:,:)].*Coeff1 + ...
    tracer(2:end,:,:).*Coeff2 + ...
    [zeros(size(tracer(1:end-1,1,:))),tracer(1:end-1,1:end-1,:)] .* Coeff3  + ...
    [tracer(1:end-1,2:end,:),zeros(size(tracer(1:end-1,1,:)))].* Coeff4;

%bottomflux=(tracer(end,:,:)-tracer(end-1,:,:))*BottomCoeff;
