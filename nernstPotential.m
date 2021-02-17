function [Erev_cell, Erev_lch, psbl_cell, psbl_lch] = nernstPotential(Cm,temp,aH2,aO2)
    global Eo R z F gamma e
    
    psbl_cell = ones(2,11);
    psbl_lch = ones(1,11);%applies only to rxns with aqueous reactants
    %Nernst potentials - electrowinning cell
    Erev_Cu_cell = Eo(1) - R*temp/(z(1)*F)*log(1/(gamma(1)*max(Cm(1),e)));
    Erev_Sn_cell = Eo(2) - R*temp/(z(2)*F)*log(1/(gamma(2)*max(Cm(2),e)));
    Erev_Al_cell = Eo(3) - R*temp/(z(3)*F)*log(1/(gamma(3)*max(Cm(3),e)));
    Erev_Pb_cell = Eo(4) - R*temp/(z(4)*F)*log(1/(gamma(4)*max(Cm(4),e)));
    Erev_Fe1_cell = Eo(5) - R*temp/(z(5)*F)*log(gamma(5)*max(Cm(5),e)/(gamma(6)*max(Cm(6),e)));
    Erev_Fe2_cell = Eo(6) - R*temp/(z(6)*F)*log(1/(gamma(5)*max(Cm(5),e)));
    Erev_Ag_cell = Eo(7) - R*temp/(z(7)*F)*log((gamma(12)*max(Cm(12),eps))^2/(gamma(7)*max(Cm(7),e)));
    Erev_Au_cell = Eo(8) - R*temp/(z(8)*F)*log((gamma(12)*max(Cm(12),eps))^2/(gamma(8)*max(Cm(8),e)));
    Erev_Pd_cell = Eo(9) - R*temp/(z(9)*F)*log((gamma(12)*max(Cm(12),eps))^4/(gamma(9)*max(Cm(9),e)));
    Erev_H_cell = Eo(10) - R*temp/(z(10)*F)*log(aH2/(gamma(10)*max(Cm(10),e))^2);
    Erev_An_cell = Eo(11) - R*temp/(z(11)*F)*log(1/(aO2*(gamma(10)*max(Cm(10),e))^4));
    Erev_cell = [Erev_Cu_cell Erev_Sn_cell Erev_Al_cell Erev_Pb_cell Erev_Fe1_cell...
    Erev_Fe2_cell Erev_Ag_cell Erev_Au_cell Erev_Pd_cell Erev_H_cell Erev_An_cell];

    if max(Cm(1),e) == e
        psbl_cell(1,1) = 0;
    end
    if max(Cm(2),e) == e
        psbl_cell(1,2) = 0;
    end
    if max(Cm(3),e) == e
        psbl_cell(1,3) = 0;
        end
    if max(Cm(4),e) == e
        psbl_cell(1,4) = 0;
    end
    if max(Cm(6),e) == e
        psbl_cell(1,6) = 0;
    end
    if max(Cm(7),e) == e
        psbl_cell(1,7) = 0;
    end
    if max(Cm(8),e) == e
        psbl_cell(1,8) = 0;
    end
    if max(Cm(9),e) == e
        psbl_cell(1,9) = 0;
    end
    if max(Cm(34),e) == e
        psbl_cell(2,1) = 0;
    end
    if max(Cm(35),e) == e
        psbl_cell(2,2) = 0;
    end
    if max(Cm(36),e) == e
        psbl_cell(2,3) = 0;
    end
    if max(Cm(37),e) == e
        psbl_cell(2,4) = 0;
    end
    if max(Cm(38),e) == e
        psbl_cell(2,6) = 0;
    end
    if max(Cm(39),e) == e
        psbl_cell(2,7) = 0;
    end
    if max(Cm(40),e) == e
        psbl_cell(2,8) = 0;
    end
    if max(Cm(41),e) == e
        psbl_cell(2,9) = 0;
    end
    
    %Nernst potentials - leaching vessel
    Erev_Cu_lch = Eo(1) - R*temp/(z(1)*F)*log(1/(gamma(1)*max(Cm(13),e)));
    Erev_Sn_lch = Eo(2) - R*temp/(z(2)*F)*log(1/(gamma(2)*max(Cm(14),e)));
    Erev_Al_lch = Eo(3) - R*temp/(z(3)*F)*log(1/(gamma(3)*max(Cm(15),e)));
    Erev_Pb_lch = Eo(4) - R*temp/(z(4)*F)*log(1/(gamma(4)*max(Cm(16),e)));
    Erev_Fe1_lch = Eo(5) - R*temp/(z(5)*F)*log(gamma(5)*max(Cm(17),e)/gamma(6)*max(Cm(18),e));
    Erev_Fe2_lch = Eo(6) - R*temp/(z(6)*F)*log(1/(gamma(5)*max(Cm(17),e)));
    Erev_Ag_lch = Eo(7) - R*temp/(z(7)*F)*log((gamma(12)*max(Cm(24),eps))^2/(gamma(7)*max(Cm(19),e)));
    Erev_Au_lch = Eo(8) - R*temp/(z(8)*F)*log((gamma(12)*max(Cm(24),eps))^2/(gamma(8)*max(Cm(20),e)));
    Erev_Pd_lch = Eo(9) - R*temp/(z(9)*F)*log((gamma(12)*max(Cm(24),eps))^4/(gamma(9)*max(Cm(21),e)));
    Erev_H_lch = Eo(10) - R*temp/(z(10)*F)*log(aH2/(gamma(10)*max(Cm(22),e))^2);
    Erev_An_lch = Eo(11) - R*temp/(z(11)*F)*log(1/(aO2*(gamma(10)*max(Cm(22),e))^4));
    Erev_lch = [Erev_Cu_lch Erev_Sn_lch Erev_Al_lch Erev_Pb_lch Erev_Fe1_lch...
        Erev_Fe2_lch Erev_Ag_lch Erev_Au_lch Erev_Pd_lch Erev_H_lch Erev_An_lch];
    
    if max(Cm(13),e) == e
        psbl_lch(1) = 0;
    end
    if max(Cm(14),e) == e
        psbl_lch(2) = 0;
    end
    if max(Cm(15),e) == e
        psbl_lch(3) = 0;
        end
    if max(Cm(16),e) == e
        psbl_lch(4) = 0;
    end
    if max(Cm(18),e) == e
        psbl_lch(6) = 0;
    end
    if max(Cm(19),e) == e
        psbl_lch(7) = 0;
    end
    if max(Cm(20),e) == e
        psbl_lch(8) = 0;
    end
    if max(Cm(21),e) == e
        psbl_lch(9) = 0;
    end
end