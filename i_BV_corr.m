%BV equation
%overpot = Electrode overpotential (V)
%i0 = Exchange current density (Use same units as desired for output current density i)
%a = charge transfer coefficient
%z = # electrons transferred
%temp = temperature (Kelvin)
function [i] = i_BVcorr(Ecorr, Erev, i0, alpha, z, temp)
	R = 8.314; %J/mol/K
    F = 96485.3329; %C/mol e-
    i = i0*(exp((alpha*z*F)/(R*temp)*(Ecorr-Erev))-exp(-((1-alpha)*z*F)/(R*temp)*(Ecorr-Erev)));
end
