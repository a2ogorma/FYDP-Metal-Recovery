yy = [resultsEconomicFeasible.metrics];
y = [yy.netAnnualafterTax];
xb = [paramSetBaseFeasible.tfinal];
xp = [paramSetPreciousFeasible.tfinal];
for k = 1:1:size(y,2)
    if xb(k) > 200000
        xbind(k) = 1;
    else
        xbind(k) = 0;
    end
    if xp(k) > 200000
        xpind(k) = 1;
    else
        xpind(k) = 0;
    end
end
xb3tot = sum(xbind);
xb2tot = size(xbind,2) - xb3tot;
xp3tot = sum(xpind);
xp2tot = size(xpind,2) - xp3tot;
xb3p3tot = sum(xbind.*xpind);
xb2p3tot = sum((-1.*(xbind-1).*xpind));
xb3p2tot = sum((-1.*(xpind-1).*xbind));
xb2p2tot = sum((-1.*(xbind-1)).*(-1.*(xpind-1)));

valb3day = sum(xbind.*y)/xb3tot;
valb2day = sum(-1.*(xbind-1).*y)./xb2tot;
valp3day = sum(xpind.*y)./xp3tot;
valp2day = sum(-1.*(xpind-1).*y)./xp2tot;
valb3p3day = sum(xbind.*xpind.*y)/xb3p3tot;
valb2p3day = sum(-1.*(xbind-1).*xpind.*y)/xb2p3tot;
valb3p2day = sum(-1.*(xpind-1).*xbind.*y)/xb3p2tot;
valb2p2day = sum((-1.*(xbind-1)).*(-1.*(xpind-1)).*y)/xb2p2tot;