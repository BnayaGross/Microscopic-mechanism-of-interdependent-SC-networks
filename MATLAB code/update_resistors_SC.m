function [R,Rin,states,states_in] = update_resistors_SC(NR, W, R, Rin, heat, heat_in, Ic, Ic_in, Tc, Tcin, Rhot, RSc, currents,cin,T0,states,states_in)

    for i=1:NR
        if (T0+heat(i))>Tc(i)
            R(i) = Rhot;
            states(i) = 3;
        else
            Icc = Ic(i)*(1-(T0+heat(i))/Tc(i))^2;
            Vres = currents(i)*R(i);    
            if Vres > Icc*Rhot
                  R(i) = Rhot;
                  states(i) = 3;
            elseif Vres < Icc*RSc
                  R(i) = RSc;
                  states(i) = 1;
            else
                R(i) = Vres/Icc;
                states(i) = 2;
            end
        end
    end

    for i=1:W
        if (T0+ heat_in(i))>Tcin(i)
            Rin(i) = Rhot;
            states_in(i) = 3;
        else
            Icc = Ic_in(i)*(1-(T0+heat_in(i))/Tcin(i))^2;
            Vres = cin(i)*Rin(i);
            if Vres > Icc*Rhot
                Rin(i) = Rhot;
                states_in(i) = 3;
            elseif Vres < Icc*RSc
                Rin(i) = RSc;
                states_in(i) = 1;
            else
                Rin(i) = Vres/Icc;  
                states_in(i) = 2;
            end
        end
    end

end
