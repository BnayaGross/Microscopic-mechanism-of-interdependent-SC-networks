% This script creates two interdependent superconducting networks, keeping
% them at a constant temperature and decrese the biased current. An aburpt
% transition and plateau behaviour is observed

% general structure of the networks
W = 100; % width size
L = 100; % length size
NR = (2*W-1)*L; % number of links
NN = W*L; % number of nodes
epsilon = 10^(-5); % threshold for epsilon to determine system converges
kappa = 10*10^6; % kappa for the media between layers
kappa_own = 5*10^1; % kappa for the media in each layer
Rhot = 500; % the resistance in the Normal phase
RSc = 10^(-5); % the resistance in the SC phase
Ic0_1 = 58*10^(-6); % avg critical current layer1
Ic0_2 = 48*10^(-6); % avg critical current layer2
sigma = 0.1; % Level of disorder
Dt = L^2; % Dt sets the level of heat diffusion. for meanfiled Dt~r^2
T = 1; % permanent temperature
Ic = nan; %initializing the critical current

% currents at entrence node properties
Ii = 200*10^(-6);
If = 500*10^(-6);
Ijumps = 1000;
I_range = linspace(Ii,If,Ijumps);
I_range = flip(I_range); % decreasing current range
LI = length(I_range);
R1_all = zeros(1,LI);
R2_all = zeros(1,LI);

%seting the general currents vector
I1 = zeros(1,NN+1);
I2 = zeros(1,NN+1);

% vector currents for each link
cin1 = zeros(1,W);
currents1 = zeros(1,NR);
cin2 = zeros(1,W);
currents2 = zeros(1,NR);

% Ic for each resistor
Ic_in1 = Ic0_1*(1+sigma*randn(1,W));
Ic1 = Ic0_1*(1+sigma*randn(1,NR));
Ic_in2 = Ic0_2*(1+sigma*randn(1,W));
Ic2 = Ic0_2*(1+sigma*randn(1,NR));

% Tc for each resistor
factor = 100;  % the 100 here can be modified but should be close to it to match the experiment
Tc1 = factor*Rhot*Ic1;
Tcin1 = factor*Rhot*Ic_in1;
Tc2 = factor*Rhot*Ic2;
Tcin2 = factor*Rhot*Ic_in2;

% vector resistors for each link. we start with all at the Normal phase
Rin1 = Rhot*ones(1,W);
R1 = Rhot*ones(1,NR);
Rin2 = Rhot*ones(1,W);
R2 = Rhot*ones(1,NR);

% vector states for each link. 1=SC , 2=intermidiate , 3=normal.
states_in1 = zeros(1,W);
states1 = zeros(1,NR);
states_in2 = zeros(1,W);
states2 = zeros(1,NR);

for i=1:length(I_range)

    I = I_range(i); % Totall currnet I
    I1(end) = I;
    I2(end) = I;

    iter = 0;
    flag = 0;                
    iter_per_I = [];
    R1_per_I = [];
    R2_per_I = [];

    while flag==0 
        iter = iter+1;

        %saving 2 last voltege results for exit condition
        if iter>2
            V1old2 = V1old;
            V2old2 = V2old;
        end

        if iter>1
            V1old = V1;
            V2old = V2;
        end

        % first, calculate the voltage at each node in the two networks    
        V1 = Voltage_nodes(L, W, NN, NR,Rin1, R1,I1);
        V2 = Voltage_nodes(L, W, NN, NR,Rin2, R2,I2);

        % checking exit condition
        if iter>2
            if norm(V1-V1old2)/norm(V1)+norm(V2-V2old2)/norm(V2)<epsilon || norm(V1-V1old)/norm(V1)+norm(V2-V2old)/norm(V2)<epsilon
                flag=1;
            end
        end

        % now, calculate the currents in each resistor and its heat production
        [currents1,cin1,heat1,heat_in_1] = find_currents(V1, L, W, NN, NR, Rin1, R1,kappa_own);
        [currents2,cin2,heat2,heat_in_2] = find_currents(V2, L, W, NN, NR, Rin2, R2,kappa_own);                

        % let the heat diffuse (only in the bulk, no exit/entrence nodes)       
        Q12 = heat_transfer_spectral(L,W,Dt,heat1(1:NR-W)*kappa);
        Q21 = heat_transfer_spectral(L,W,Dt,heat2(1:NR-W)*kappa);
        Q11 = heat_transfer_spectral(L,W,100,heat1(1:NR-W));
        Q22 = heat_transfer_spectral(L,W,100,heat2(1:NR-W));

        heat1(1:NR-W) = Q21+Q11;                  
        heat2(1:NR-W) =  Q12+Q22;

        [R1 , Rin1, states1 ,states_in1] = update_resistors_SC(NR, W, R1, Rin1, heat1, heat_in_1, Ic1, Ic_in1, Tc1, Tcin1, Rhot, RSc, currents1, cin1,T, states1 ,states_in1);
        [R2 , Rin2, states2 ,states_in2] = update_resistors_SC(NR, W, R2, Rin2, heat2, heat_in_2, Ic2, Ic_in2, Tc2, Tcin2, Rhot, RSc, currents2, cin2,T, states2 ,states_in2);

        R1_per_I = [R1_per_I V1(end)/I];
        R2_per_I = [R2_per_I V2(end)/I];
        iter_per_I = [iter_per_I iter];
        
    end
    
    % updating resistence
    R1_all(i) = V1(end)/I;
    R2_all(i) = V2(end)/I;
    
    % keep the process data if we are at the critical first order phase transition point 
    if abs(R1_per_I(1)-R1_per_I(end)) > 200 && abs(R2_per_I(1)-R2_per_I(end)) > 200
        R1_plateau = R1_per_I;
        R2_plateau = R2_per_I;
        iteration_plateau = iter_per_I;  
        Ic = I;
    end
    
end

%resistance plot
figure;
plot(I_range, R1_all,'-r*', I_range, R2_all,'-b*');
title('Resistence by lower the current');
legend('layer1','layer2', 'Location', 'best');
ylabel('R(\Omega)');
xlabel('I(A)'); 

%plateau plot
if ~isnan(Ic)
    figure;
    plot(iteration_plateau, R1_plateau,'-r*', iteration_plateau, R2_plateau,'-b*');
    title(['Resistence at criticality, ', 'I_c = ', num2str(Ic*10^(6)), ' [\muA]']);
    legend('layer1','layer2', 'Location', 'best');
    ylabel('R(\Omega)');
    xlabel('iteration'); 
end

