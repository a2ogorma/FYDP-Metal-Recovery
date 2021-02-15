%BV equation
%overpot = Electrode overpotential (V)
%i0 = Exchange current density (Use same units as desired for output current density i)
%a = charge transfer coefficient
%z = # electrons transferred
%temp = temperature (Kelvin)d
%add situation 1 - overpotential 2 3
function i_f = i_BV(overpot, i0, iL, a, z, temp)
    R = 8.314; %J/mol/K
    F = 96485.3329; %C/mol e-
    j_BV = zeros(1,numel(overpot));
    for idx = 1:1:numel(overpot)
        if overpot(idx)<0.15 && overpot(idx)>-0.15
            %Total kinetic control, no Tafel approximation
            j_BV(idx) = i0(idx)*(exp((a(idx)*z(idx)*F)/(R*temp)*(overpot(idx)))...
                -exp(-((1-a(idx))*z(idx)*F)/(R*temp)*(overpot(idx))));
        elseif overpot(idx)>=0.15
            %Anodic mass-transfer kinetic mixed case, Tafel approx
            rhs = exp(log(i0(idx)/iL(idx)/z(idx))+a(idx)*z(idx)*F*overpot(idx)/(R*temp));
            j_BV(idx) = iL(idx)*rhs/(1+rhs);
        elseif overpot(idx)<=-0.15
            %cathodic mass-transfer kinetic mixed case, Tafel approx
            rhs = exp(log(i0(idx)/iL(idx)/z(idx))-(1-a(idx))*z(idx)*F*overpot(idx)/(R*temp));
            j_BV(idx) = -(iL(idx)*rhs)/(1+rhs);
        end
    end
    i_f = j_BV;
end