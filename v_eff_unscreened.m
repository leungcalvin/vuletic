function v_eff_vals = v_eff_unscreened(r,l_oam,netZ);
    % Generates U_eff in either S.I. or natural units where V(r) = Ze^2/(4*pi*eps0*r)
    RBOHR=  5.29177210*10^(-11); %m
    alpha = 1/137;
    HBAR = 1.054e-27; %erg-s
    MASS =9.10938215e-28; %electron mass in grams
    %V = netZ * HBAR.^2 ./ (MASS * RBOHR .* r) ;
    %L = @(l) l*(l+1)*HBAR.^2 ./(2 .* MASS .* r.^2);
    %DEMO ARGS
    eps = RBOHR * 0.01;
    r = (eps:eps:5*RBOHR)';
    netZ = 2;
    Veff_si_units =  @(l) l*(l+1)*HBAR.^2 ./(2 .* MASS .* r.^2) - netZ * HBAR.^2 ./ (MASS * RBOHR .* r);
    Veff_mikami_units = @(l) l*(l+1)* alpha * RBOHR ./(2 .* r.^2)  - netZ * alpha ./ r;
    %DEMO PLOT

    plot(r/RBOHR,0*r,'black');hold on;
    plot(r/RBOHR,Veff_si_units(0),'DisplayName','veff0');
    plot(r/RBOHR,Veff_si_units(1),'DisplayName','veff1');
    plot(r/RBOHR,Veff_si_units(2),'DisplayName','veff2');
    plot(r/RBOHR,Veff_si_units(3),'DisplayName','veff3');
    legend('-DynamicLegend');

    v_eff_vals = Veff_si_units(l_oam);
    
end
