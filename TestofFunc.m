Vapp = 10;
clear t eta1 eta2 value
eta1(1) = 0;
eta2(1) = -0.1;
value(1) = 0;
for iter = 1:1:1000
    t(iter) = iter;
    CurrentFunc = @(eta2) Aan*( iAn0*(exp(alphaAn*4*F*eta1(iter)/(Rgas*T))-exp(-(1-alphaAn)*4*F*eta1(iter)/(Rgas*T))) + iFe0*(exp(alphaFe*F*(-ErevAn(1)+eta1(iter)+ErevFe(1))/(Rgas*T))-exp(-(1-alphaFe)*F*(-ErevAn(1)+eta1(iter)+ErevFe(1))/(Rgas*T))))+Acat*( iAg0*(exp(alphaAg*F*eta2/(Rgas*T))-exp(-(1-alphaAg)*F*eta2/(Rgas*T))) + iAu0*(exp(alphaAu*F*(ErevAg(1)+eta2-ErevAu(1))/(Rgas*T))-exp(-(1-alphaAu)*F*(ErevAg(1)+eta2-ErevAu(1))/(Rgas*T))) + iPd0*(exp(alphaPd*2*F*(ErevAg(1)+eta2-ErevPd(1))/(Rgas*T))-exp(-(1-alphaPd)*2*F*(ErevAg(1)+eta2-ErevPd(1))/(Rgas*T))));
    eta2(iter+1) = fzero(CurrentFunc,eta2(iter));
    eta1(iter+1) = eta1(iter)+0.01;
    value(iter+1) = CurrentFunc(eta2(iter));
end
plot(t,eta1,t,eta2,t,value)
