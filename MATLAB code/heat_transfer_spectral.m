function [Q] = heat_transfer_spectral(L,W,Dt,Q)

% Dt = 0.2;
% L=30;
% W=30;
% NR1 = L*(W-1);
% NR2 = W*(L-1);
% Q = zeros(1,NR1+NR2);
% Q(round((NR2+NR1)/4)+L/2) = 1;

% takes the L and the W plus a vector Q for the heat in each pair of
% resistors and build a heat matrix.

% these will be the parameters of the function, Q is the vector of
% heat for each bulk resistor

% NR1 = L*(W-1);
% NR2 = W*(L-1);
% Q = rand(1,NR1+NR2);
%%% end parameters

% build the matrix

heat = zeros(2*W-1,2*L-1);

mone = 1;
for j=1:2:2*L-1
    for i=2:2:2*W-1
        heat(i,j) = Q(mone);
        mone = mone+1;
    end
end

for j=2:2:2*L-1
    for i=1:2:2*W-1
        heat(i,j) = Q(mone);
        mone = mone+1;
    end
end


% this is the heat map before diffusion, odd sites are empty
% figure(22); imagesc(heat);
% colorbar;

% this measures the total heat before diffusion
%sum(sum(heat))

% diffusion of heat
% heat = fft2d_heat(Dt,heat);
heat = fft2d(Dt,heat);

% total heat after diffusion
%sum(sum(heat))

% this is the heat map after diffusion, odd sites are hot too
% figure; imagesc(heat);
% colorbar;

% takes the heat out of the empty sites
heat1 = zeros(2*W+1,2*L+1);

for i = 1:2:2*W-1
    for j = 1:2:2*L-1
        heat1(1+i-1,1+j) = heat1(1+i-1,1+j)+heat(i,j)/4;
        heat1(1+i+1,1+j) = heat1(1+i+1,1+j)+heat(i,j)/4;
        heat1(1+i,1+j-1) = heat1(1+i,1+j-1)+heat(i,j)/4;
        heat1(1+i,1+j+1) = heat1(1+i,1+j+1)+heat(i,j)/4;
        heat(i,j)=0;
    end
end

for i = 2:2:2*W-1
    for j = 2:2:2*L-1
        heat1(1+i-1,1+j) = heat1(1+i-1,1+j)+heat(i,j)/4;
        heat1(1+i+1,1+j) = heat1(1+i+1,1+j)+heat(i,j)/4;
        heat1(1+i,1+j-1) = heat1(1+i,1+j-1)+heat(i,j)/4;
        heat1(1+i,1+j+1) = heat1(1+i,1+j+1)+heat(i,j)/4;
        heat(i,j)=0;
    end
end

heat = heat+heat1(2:end-1,2:end-1);

for i=2:2:2*W+1
    heat(i-1,end-1)= heat(i-1,end-1)+heat1(i,1);
    heat(i-1,2)= heat(i-1,2)+heat1(i,end);
   
end

for j=2:2:2*L+1
    heat(end-1,j-1)= heat(end-1,j-1)+heat1(1,j);
    heat(2,j-1)= heat(2,j-1)+heat1(end,j);
end


% this is a heat map after the heat of the odd site has been taken to the even
%figure(24); imagesc(heat)
%colorbar

% returns a heat to every resistor

mone = 1;
for j=1:2:2*L-1
    for i=2:2:2*W-1
        Q(mone) = heat(i,j);
        mone = mone+1;
    end
end

for j=2:2:2*L-1
    for i=1:2:2*W-1
        Q(mone) = heat(i,j);
        mone = mone+1;
    end
end

% total heat in the return map  
%sum(Q)

 
 
