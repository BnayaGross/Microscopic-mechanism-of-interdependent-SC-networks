function [current,current_in,heat,heat_in] = find_currents(V, L, W, NN, NR, Rin, R,kappa)

heat = zeros(1,NR);
heat_in = zeros(1,W);
current = zeros(1,NR);
current_in = zeros(1,W);

mone = 1;
% calculate curents and dissipated heat - columns

for k=1:L
    for i=1:W-1
        m = W*(k-1)+i;
        J = (V(m)-V(m+1))/R(mone);
        heat(mone) = kappa*J^2*R(mone);
        current(mone) = J;
        mone = mone +1;
    end
end

% calculate currents and dissipated heat - rows

for k=1:L-1
    for i=1:W
        m1 = W*(k-1)+i;
        m2 = W*k+i;
        J = (V(m1)-V(m2))/R(mone);
        heat(mone) = kappa*J^2*R(mone);   
        current(mone) = J;
        mone = mone +1;
    end
end

% calculate cuurents and dissipated heat - exits
 
 for i=1:W
     J = V(NN-W+i)/R(mone);
     heat(mone) = kappa*J^2*R(mone);
     current(mone) = J;
     mone = mone+1;
 end
     
% calculate cuurents and dissipated heat - entrence 

 for i=1:W
     J = (V(end)-V(i))/Rin(i);
    current_in(i) = J;
    heat_in(i) = kappa*J^2*Rin(i);  
 end
