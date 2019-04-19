function powercheck(x,y,refpointidx,power)
    figure;
    loglog(x,y,'DisplayName','y(x)');hold on;
    powerlaw = y(refpointidx)*(x/x(refpointidx)).^power;
    powerlaw(powerlaw > max(y)) = max(y);
    powerlaw(powerlaw < min(y)) = min(y);
    plot(x,powerlaw,'DisplayName',string(power));
    legend('-DynamicLegend');
    