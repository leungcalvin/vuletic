function v_eff_vals = v_eff_unscreened2(r,l_oam,netZ);
    % Generates V_eff in units where e^2 = alpha; 4pieps0 = 1; hbar = 1;
    % V(r) = Ze^2/(4*pi*eps0*r)
    RBOHR=  5.29177210*10^(-11); %m
    alpha = 1/137;

    %DEMO
%     eps = RBOHR * 0.01;
%     r = (eps:eps:5*RBOHR)';
%     netZ = 2;
    Veff =  @(l) l*(l+1) * RBOHR * alpha ./(2.* r.^2) - netZ * alpha ./ r;
    %DEMO

%     plot(r/RBOHR,0*r,'black');hold on;axis([0 ,5 ,-2e8,3e8]);
%     plot(r/RBOHR,Veff(0),'DisplayName','veff0');
%     plot(r/RBOHR,Veff(1),'DisplayName','veff1');
%     plot(r/RBOHR,Veff(2),'DisplayName','veff2');
%     plot(r/RBOHR,Veff(3),'DisplayName','veff3');
    legend('-DynamicLegend');

    v_eff_vals = Veff(l_oam);
    
end
