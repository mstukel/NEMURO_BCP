function[profile_new,Flux,Coeff1,Coeff2]=sinking2D_ftcs(profile,Coeff1,Coeff2)

% thickness=z_int(2:end)-z_int(1:end-1);
% 
% for j=1:length(profile(1,:))
%     i=1;
%     Coeff1(i,j) = 1 - omega(j)*dt/thickness(i);
%     Coeff2(i,j) = 0;
%     Coeff3(i,j) = omega(j)*dt;
%     for i=2:length(profile(:,1))
%         Coeff1(i,j) = 1 - omega(j)*dt/thickness(i);
%         Coeff2(i,j) = omega(j)*dt/thickness(i);
%         Coeff3(i,j) = omega(j)*dt;
%     end
% end

profile_new = profile.*Coeff1 + [zeros(1, length(profile(1,:,1)), length(profile(1,1,:))); profile(1:end-1,:,:)].*Coeff2;

%Flux=profile.*Coeff3;

