yy = [resultsEconomic.metrics];
y = [yy.netAnnualafterTax];
ytitle = 'Net Annual Profit after tax (CA$)';

%{
for k = 1:10
    figure
    scatter(x(k).valsPrec,y)
    xlabel(x(k).name)
    ylabel(ytitle)
    title('Precious Metal')
end
for k = 1:10
    figure
    scatter(x(k).valsBase,y)
    xlabel(x(k).name)
    ylabel(ytitle)
    title('Base Metal')
end
%}