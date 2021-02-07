%BV equation
%overpot = Electrode overpotential (V)
%i0 = Exchange current density (Use same units as desired for output current density i)
%a = charge transfer coefficient
%z = # electrons transferred
%temp = temperature (Kelvin)
function [i_f] = i_BV(overpot, i0, a, z, temp)
	R = 8.314; %J/mol/K
    F = 96485.3329; %C/mol e-
    j_BV = i0.*(exp((a.*z*F)/(R*temp).*(overpot))-exp(-((1-a).*z*F)/(R*temp).*(overpot)));
    for k = 1:1:numel(j_BV)
        if j_BV(k)>(realmax/1E20)
            j_BV(k) = realmax/1E20;
        elseif j_BV(k)< (-realmax/1E20)
            j_BV(k) = -realmax/1E20;
        end
    end
    i_f = j_BV;
end