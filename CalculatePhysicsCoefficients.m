function[Coeff0,Coeff1,Coeff2]=CalculatePhysicsCoefficients(tracer,z,z_edge,Kz,w,dt,deep)

%This function is written to calculate the coefficients necessary for
%one-dimensional advection and diffusion.  It is NOT conservative if there
%is divergence in the advection terms.  It is paired with physics_ftcs.

%Inputs:
%tracer = tracer concentrations at grid center points (z).  Note that tracer is only used to define the size of the coefficients, which then allows more rapid computations in physics_ftcs
%z = depths of the midpoints of the grid cells
%z_edge = depths of the edges of the grid cells
%Kz = vertical diffusivity at the edges of the grid cells (should be 0 at the surface)
%w = vertical velocities at the edges of the grid cells (should be 0 at the surface), negative is treated as a downward velocity
%dt is the model's physical time step
%deep is the tracer concentrations immediately below the model (boundary conditions)

%Note that this code is written to be transparent, rather than
%computationally efficient, because it is assumed that this code is run
%once and that physics_ftcs is then run repeatedly.


deltaz=z(2:end)-z(1:end-1);
deltaz_edge=z_edge(2:end)-z_edge(1:end-1);

numtracers = size(tracer);
numtracers = numtracers(2);

tracer(end+1,:)=deep;


for j=1:numtracers
    i=1;
    MixCoeff0(i,j) = 1 - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
    MixCoeff1(i,j) = 0;
    MixCoeff2(i,j) = (Kz(i+1)*dt)./(deltaz_edge(i)*deltaz(i));
    Coeff0(i,j) = 1 - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
    for i=2:length(tracer(:,1))-1
        MixCoeff0(i,j) = 1 - (Kz(i) * dt)./(deltaz_edge(i).*deltaz(i-1)) - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
        MixCoeff1(i,j) = (Kz(i)*dt)./(deltaz_edge(i)*deltaz(i-1));
        MixCoeff2(i,j) = (Kz(i+1)*dt)./(deltaz_edge(i)*deltaz(i));
        Coeff0(i,j) = 1 - (Kz(i) * dt)./(deltaz_edge(i).*deltaz(i-1)) - (Kz(i+1) *dt)./(deltaz_edge(i)*deltaz(i));
    end
end


% tracernew = tracer(1:end-1,:).*MixCoeff0 + ...
%     [zeros(1,length(tracer(1,:)));tracer(1:end-2,:)].*MixCoeff1 + ...
%     tracer(2:end,:).*MixCoeff2;

w_pos = max(w,0);
w_neg = min(w,0);


for j=1:numtracers
    i=1;
    AdvDownCoeff0(i,j) = 1;
    AdvDownCoeff1(i,j) = 0;
    for i=2:length(tracer(:,1))-1
        AdvDownCoeff0(i,j) = 1 + w_neg(i)*dt/deltaz_edge(i);
        AdvDownCoeff1(i,j) = -w_neg(i)*dt/deltaz_edge(i);
        Coeff0(i,j) = Coeff0(i,j) + w_neg(i)*dt/deltaz_edge(i);
    end
end

% profile_new = profile.*Coeff0 + [zeros(1,length(w_neg));profile(1:end-1,:)].*Coeff1;


for j=1:numtracers
    for i=1:length(tracer(:,1))-1
        AdvUpCoeff0(i,j) = 1 - w_pos(i+1)*dt/deltaz_edge(i);
        AdvUpCoeff2(i,j) = w_pos(i+1)*dt/deltaz_edge(i);
        Coeff0(i,j) = Coeff0(i,j) - w_pos(i+1)*dt/deltaz_edge(i);
    end
end

Coeff1 = MixCoeff1 + AdvDownCoeff1;
Coeff2 = MixCoeff2 + AdvUpCoeff2;


% profile_new = profile.*Coeff0 + [profile(2:end,:);Profile(end,:)].*Coeff2;