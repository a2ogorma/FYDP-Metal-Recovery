function [Erev_cat, Erev_an, Erev_lch] = nernstPotential(Cm,temp,solution)
    global Eo R z F gamma e aH2 aO2
    %{
    Reactions
    Cu2+ + 2e- <--> Cu(s) (1)
    Sn2+ + 2e- <--> Sn(s) (2)
    Fe3+ + e- <--> Fe2+ (3)
    Fe2+ + 2e- <--> Fe(s) (4)
    Ag+ + e- <--> Ag(s) (5B) OR [Ag(S2O3)2]3- + e- <--> Ag(s) + 2(S2O3)2- (5P)
    AgCl(s) + e- <--> Ag(s) + Cl- (6)
    Au3+ + 3e- <--> Au(s) (7B) OR [Au(S2O3)2]3- + e- <--> Au(s) + 2(S2O3)2- (7P)
    [AuCl4]- + 3e- <--> Au(s) + 4Cl- (8)
    Pd2+ + 2e- <--> Pd(s) (9B) OR [Pd(S2O3)4]6- + 2e- <--> Pd(s) + 4(S2O3)2- (9P)
    2H+ + e- <--> H2(g) (10)
    4H+ + O2(g) + 4e- <--> 2H2O(l) (11)
    %}
    
    psbl_cat = ones(1,11);
    psbl_an = ones(1,11);
    psbl_lch = ones(1,11);%applies only to rxns with aqueous reactants
    
    %Nernst potentials - electrowinning cell
    %Cathode:
    Erev_1_cat = Eo(1) - R*temp/(z(1)*F)*log(1/(gamma(1)*max(Cm(1),e)));
    Erev_2_cat = Eo(2) - R*temp/(z(2)*F)*log(1/(gamma(2)*max(Cm(2),e)));
    Erev_3_cat = Eo(3) - R*temp/(z(3)*F)*log((gamma(3)*max(Cm(3),e))/(gamma(4)*max(Cm(4),e)));
    Erev_4_cat = Eo(4) - R*temp/(z(4)*F)*log(1/(gamma(3)*max(Cm(3),e)));
    Erev_6_cat = Eo(6) - R*temp/(z(6)*F)*log(gamma(9)*max(Cm(9),e));
    Erev_8_cat = Eo(8) - R*temp/(z(8)*F)*log((gamma(9)*max(Cm(9),e))^4/(gamma(10)*max(Cm(10),e)));
    Erev_10_cat = Eo(10) - R*temp/(z(10)*F)*log(aH2^0.5/(gamma(8)*max(Cm(8),e)));
    Erev_11_cat = Eo(11) - R*temp/(z(11)*F)*log(1/(aO2^0.25*(gamma(8)*max(Cm(8),e))));
    %Variable reactions:
    if solution == 1 %Base metal system
        Erev_5_cat = Eo(5) - R*temp/(z(5)*F)*log(1/(gamma(5)*max(Cm(5),e)));
        Erev_7_cat = Eo(7) - R*temp/(z(7)*F)*log(1/(gamma(6)*max(Cm(6),e)));
        Erev_9_cat = Eo(9) - R*temp/(z(9)*F)*log(1/(gamma(7)*max(Cm(7),e)));
    else %Precious metal system
        Erev_5_cat = Eo(5) - R*temp/(z(5)*F)*log((gamma(9)*max(Cm(9),e))^2/(gamma(5)*max(Cm(5),e)));
        Erev_7_cat = Eo(7) - R*temp/(z(7)*F)*log((gamma(9)*max(Cm(9),e))^2/(gamma(6)*max(Cm(6),e)));
        Erev_9_cat = Eo(9) - R*temp/(z(9)*F)*log((gamma(9)*max(Cm(9),e))^4/(gamma(7)*max(Cm(7),e)));
    end
    Erev_cat = [Erev_1_cat Erev_2_cat Erev_3_cat Erev_4_cat Erev_5_cat...
        Erev_6_cat Erev_7_cat Erev_8_cat Erev_9_cat Erev_10_cat Erev_11_cat];
    
    %Anode:
    Erev_1_an = Eo(1) - R*temp/(z(1)*F)*log(1/(gamma(1)*max(Cm(11),e)));
    Erev_2_an = Eo(2) - R*temp/(z(2)*F)*log(1/(gamma(2)*max(Cm(12),e)));
    Erev_3_an = Eo(3) - R*temp/(z(3)*F)*log((gamma(3)*max(Cm(13),e))/(gamma(4)*max(Cm(14),e)));
    Erev_4_an = Eo(4) - R*temp/(z(4)*F)*log(1/(gamma(3)*max(Cm(13),e)));
    Erev_6_an = Eo(6) - R*temp/(z(6)*F)*log(gamma(9)*max(Cm(19),e));
    Erev_8_an = Eo(8) - R*temp/(z(8)*F)*log((gamma(9)*max(Cm(19),e))^4/(gamma(10)*max(Cm(20),e)));
    Erev_10_an = Eo(10) - R*temp/(z(10)*F)*log(aH2^0.5/(gamma(8)*max(Cm(18),e)));
    Erev_11_an = Eo(11) - R*temp/(z(11)*F)*log(1/(aO2^0.25*(gamma(8)*max(Cm(18),e))));
    %Variable reactions:
    if solution == 1 %Base metal system
        Erev_5_an = Eo(5) - R*temp/(z(5)*F)*log(1/(gamma(5)*max(Cm(15),e)));
        Erev_7_an = Eo(7) - R*temp/(z(7)*F)*log(1/(gamma(6)*max(Cm(16),e)));
        Erev_9_an = Eo(9) - R*temp/(z(9)*F)*log(1/(gamma(7)*max(Cm(17),e)));
    else %Precious metal system
        Erev_5_an = Eo(5) - R*temp/(z(5)*F)*log((gamma(9)*max(Cm(19),e))^2/(gamma(5)*max(Cm(15),e)));
        Erev_7_an = Eo(7) - R*temp/(z(7)*F)*log((gamma(9)*max(Cm(19),e))^2/(gamma(6)*max(Cm(16),e)));
        Erev_9_an = Eo(9) - R*temp/(z(9)*F)*log((gamma(9)*max(Cm(19),e))^4/(gamma(7)*max(Cm(17),e)));
    end
    Erev_an = [Erev_1_an Erev_2_an Erev_3_an Erev_4_an Erev_5_an...
        Erev_6_an Erev_7_an Erev_8_an Erev_9_an Erev_10_an Erev_11_an];
    %{
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
    %}
    %Nernst potentials - leaching vessel
    Erev_1_lch = Eo(1) - R*temp/(z(1)*F)*log(1/(gamma(1)*max(Cm(21),e)));
    Erev_2_lch = Eo(2) - R*temp/(z(2)*F)*log(1/(gamma(2)*max(Cm(22),e)));
    Erev_3_lch = Eo(3) - R*temp/(z(3)*F)*log((gamma(3)*max(Cm(23),e))/(gamma(4)*max(Cm(24),e)));
    Erev_4_lch = Eo(4) - R*temp/(z(4)*F)*log(1/(gamma(3)*max(Cm(23),e)));
    Erev_6_lch = Eo(6) - R*temp/(z(6)*F)*log(gamma(9)*max(Cm(29),e));
    Erev_8_lch = Eo(8) - R*temp/(z(8)*F)*log((gamma(9)*max(Cm(29),e))^4/(gamma(10)*max(Cm(30),e)));
    Erev_10_lch = Eo(10) - R*temp/(z(10)*F)*log(aH2^0.5/(gamma(8)*max(Cm(28),e)));
    Erev_11_lch = Eo(11) - R*temp/(z(11)*F)*log(1/(aO2^0.25*(gamma(8)*max(Cm(28),e))));
    %Variable reactions:
    if solution == 1 %Base metal system
        Erev_5_lch = Eo(5) - R*temp/(z(5)*F)*log(1/(gamma(5)*max(Cm(25),e)));
        Erev_7_lch = Eo(7) - R*temp/(z(7)*F)*log(1/(gamma(6)*max(Cm(26),e)));
        Erev_9_lch = Eo(9) - R*temp/(z(9)*F)*log(1/(gamma(7)*max(Cm(27),e)));
    else %Precious metal system
        Erev_5_lch = Eo(5) - R*temp/(z(5)*F)*log((gamma(9)*max(Cm(29),e))^2/(gamma(5)*max(Cm(25),e)));
        Erev_7_lch = Eo(7) - R*temp/(z(7)*F)*log((gamma(9)*max(Cm(29),e))^2/(gamma(6)*max(Cm(26),e)));
        Erev_9_lch = Eo(9) - R*temp/(z(9)*F)*log((gamma(9)*max(Cm(29),e))^4/(gamma(7)*max(Cm(27),e)));
    end
    Erev_lch = [Erev_1_lch Erev_2_lch Erev_3_lch Erev_4_lch Erev_5_lch...
        Erev_6_lch Erev_7_lch Erev_8_lch Erev_9_lch Erev_10_lch Erev_11_lch];
    %{
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
        %}