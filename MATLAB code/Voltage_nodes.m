function [V] = Voltage_nodes(L, W, NN, NR, Rin, R,I)
    
%%%%%%%%%% build the matrix

isparse = zeros(1,NR+W);
jsparse = zeros(1,NR+W);
Gsparse = zeros(1,NR+W);
mone = 1;

for i=1:L
    for j=1:W
        m1 = W*(i-1)+j; % node number
        
        a = (W-1)*(i-1)+j-1; % resistor number above vertex, relevant if j>1
        if j>1
            Ga = 1/R(a);
        else
            Ga = 0;
        end
        
        b = (W-1)*(i-1)+j; % resistor number below vertex,relevant if j<W
        if j<W
            Gb = 1/R(b);
        else
            Gb = 0;
        end
        
        c = L*(W-1)+W*(i-1)+j; % resistor number to the right vertex, include exit
        Gr = 1/R(c);
       
        d = L*(W-1)+W*(i-2)+j; % resistor number to the left vertex + entrance
        if i>1
            Gl = 1/R(d);
        else
            Gl = 1/Rin(j);
        end
        
        % fill diagonal
        isparse(mone) = m1;
        jsparse(mone) = m1;
        Gsparse(mone) = Ga+Gb+Gr+Gl;

        mone = mone+1;
        
        % add below
        if j<W
            m2 = m1+1;
            
            isparse(mone) = m1;
            jsparse(mone) = m2;
            Gsparse(mone) = -Gb;
            mone = mone+1;
            
            isparse(mone) = m2;
            jsparse(mone) = m1;
            Gsparse(mone) = -Gb;
            mone=mone+1;
        end

        % add to the right
        if i<L
            m2 = m1+W;
            
            isparse(mone) = m1;
            jsparse(mone) = m2;
            Gsparse(mone) = -Gr;
            mone = mone+1;
            
            isparse(mone) = m2;
            jsparse(mone) = m1;
            Gsparse(mone) = -Gr;
            mone=mone+1;
        end
    end
end

% add lower line and right column for the zero node
% first the diagonal
isparse(mone) = NN+1;
jsparse(mone) = NN+1;
Gsparse(mone) = sum(1./Rin);
mone = mone+1;

% then the off diagonals
for i=1:W
    isparse(mone) = i;
    jsparse(mone) = NN+1;
    Gsparse(mone) = -1/Rin(i);
    mone = mone+1;
    
    isparse(mone) = NN+1;
    jsparse(mone) = i;
    Gsparse(mone) = -1/Rin(i);
    mone = mone+1;
end

G = sparse(isparse,jsparse,Gsparse);

%%%%%%%%%%%%%%%%%%%%%%%

% find nodes voltage
% note that the LAST term is the voltage on the zero (entrence) node

 V = mldivide(G,I');

end