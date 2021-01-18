function [func] = cor(Ecorr, Erev, i0, alpha, z, temp)
  %1 = Cu, 2= Fe  
  i_Cu = i_BV(Ecorr-Erev(1), i0(1), alpha(1), z(1), temp);
  i_Fe = i_BV(Ecorr-Erev(2), i0(2), alpha(2), z(2), temp);
  func = i_Cu + i_Fe;
  end
