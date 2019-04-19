%%% Compute Radial Density
%l = 0, max n = 6
[r,dw,rR_0,rR_0a] = joon_tf_schrodinger(0,6,0,2000,10000);
[r,dw,rR_1,rR_1a] = joon_tf_schrodinger(1,6,0,2000,10000);
[r,dw,rR_2,rR_2a] = joon_tf_schrodinger(2,6,0,2000,10000);
[r,dw,rR_3,rR_3a] = joon_tf_schrodinger(3,6,0,2000,10000);
Rfun = @(Z,n,l,r) sqrt(Z^3*(2/n^2)^2*factorial(n-l-1)/factorial(n+l))*exp(-Z*r/n).*(2*Z*r/n).^l.*laguerreL(n-l-1,2*l+1,2*Z*r/n);

density0 = rR_0.^2;
density1 = rR_1.^2;
density2 = rR_2.^2;
density3 = rR_3.^2;

density0a = rR_0a.^2;
density1a = rR_1a.^2;
density2a = rR_2a.^2;
density3a = rR_3a.^2;


%% 935 transition (no spin coupling) (numerical+screened vs analytic)
% n = 4, l = 3, rn+1 = 1 -> n = 6, l = 0, rn+1 = 6
drRsquared935 = density3(:,1) - density0(:,6);
drRsquared935a = ((r.*Rfun(70,4,3,r)).^2 - (r.*Rfun(70,6,0,r)).^2)';

%% 399/369 transition (no spin coupling) (numerical+screened vs analytic)
% n = 6, l = 1, rn+1 = 5 -> n = 6, l = 0, rn+1 = 6
drRsquared399 = density1(:,5) - density0(:,6);
drRsquared399a = ((r.*Rfun(70,6,1,r)).^2 - (r.*Rfun(70,6,0,r)).^2)';

%% 411/435 transition (no spin coupling) (numerical+screened vs analytic)
% n = 6, l = 2, rn+1 = 4 -> n = 6, l = 0, rn+1 = 6
drRsquared411 = density2(:,4) - density0(:,6);
drRsquared411a = ((r.*Rfun(70,6,2,r)).^2 - (r.*Rfun(70,6,0,r)).^2)';


masses = logspace(1,log10(3e7),1000)'; %10 eV-1e7 eV

ps935test = particleshift2(r,dw,drRsquared935',masses);
ps935testa = particleshift2(r,dw,drRsquared935a',masses);
ps369test = particleshift2(r,dw,drRsquared399',masses);
ps369testa = particleshift2(r,dw,drRsquared399a',masses);
ps411test = particleshift2(r,dw,drRsquared411',masses);
ps411testa = particleshift2(r,dw,drRsquared411a',masses);


%ps935 = particleshift(r,dw,density3(:,1)',density0(:,6)',masses);
%ps399 = particleshift(r,dw,density1(:,5)',density0(:,6)',masses);

%% 935/369 comparison
kingslope = 2.4928;figure;power = 2;
constraint = 1./(abs(kingslope*ps935test - ps369test));
loglog(masses,constraint,'DisplayName','screened non-rel. WF');hold on;
title(['kingslope = ', string(kingslope)]);
axis([min(masses),max(masses),1e-18,1e-4]);
loglog(masses,1./(abs(kingslope*ps935testa - ps369testa)),'linestyle','--','DisplayName','analytical non-rel. WF')
loglog(masses,masses.^power .* constraint(numel(constraint))./masses(numel(masses)).^power,'DisplayName',['powerlaw exp='+string(power)])
legend('-DynamicLegend');
savefig('~/Documents/vuletic/2011228_constraint_correct.fig')

kingslope = 1;figure;power = 4;
constraint = 1./(abs(kingslope*ps935test - ps369test));
loglog(masses,constraint,'DisplayName','screened non-rel. WF');hold on;
title(['kingslope = ', string(kingslope)]);
axis([min(masses),max(masses),1e-18,1e-4]);
loglog(masses,1./(abs(kingslope*ps935testa - ps369testa)),'linestyle','--','DisplayName','analytical non-rel. WF')
loglog(masses,masses.^power .* constraint(numel(constraint))./masses(numel(masses)).^power,'DisplayName',['powerlaw exp='+string(power)])
legend('-DynamicLegend');
savefig('~/Documents/vuletic/2011228_constraint_mikami.fig')

%% 369/411 comparison
kingslope = -0.74957;figure;power = 2;
constraint = 1./(abs(kingslope*ps411test - ps369test));
loglog(masses,constraint,'DisplayName','screened non-rel. WF (1 Hz)');hold on;
loglog(masses,constraint * 9e6,'DisplayName','screened non-rel. WF (9 MHz)');hold on;

title(['kingslope = ', string(kingslope)]);
axis([min(masses),max(masses),1e-18,1e-6]);
loglog(masses,1./(abs(kingslope*ps411testa - ps369testa)),'linestyle','--','DisplayName','analytical non-rel. WF (1 Hz)')
%loglog(masses,masses.^power .* constraint(numel(constraint))./masses(numel(masses)).^power,'DisplayName',['powerlaw exp='+string(power)])
legend('-DynamicLegend');

kingslope = 1;figure;power = 4;
constraint = 1./(abs(kingslope*ps411test - ps369test));
loglog(masses,constraint,'DisplayName','screened non-rel. WF');hold on;
title(['kingslope = ', string(kingslope)]);
axis([min(masses),max(masses),1e-18,1e-6]);
loglog(masses,1./(abs(kingslope*ps411testa - ps369testa)),'linestyle','--','DisplayName','analytical non-rel. WF')
loglog(masses,masses.^power .* constraint(numel(constraint))./masses(numel(masses)).^power,'DisplayName',['powerlaw exp='+string(power)])
legend('-DynamicLegend');

