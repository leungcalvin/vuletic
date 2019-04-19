SL = 3e8; %m/s
HBAR = 1.054e-34; %J-s
RBOHR=  5.29177210*10^(-11); %m
EVperJ = 1.602*10^(-19); %1 eV / 1 Joule
Y = @(m,r) (2.*pi).^-1.*SL.*exp(-m.*EVperJ.*r ./ (HBAR* SL)) ./ (4*pi.*r);
Rhyd = @(r,n,l) sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2/(n.*RBOHR))^3)).*((2.*r)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*r)./(n *RBOHR)) .*exp(-r./(n*RBOHR));
r = linspace(0,RBOHR*50,10000);
dr = (r(2)-r(1)).*ones(size(r));
%semilogy(r,Rhyd(r,4,3));


integrand935 = @(r) r.^2 .* Y(10,r) .* (Rhyd(r,4,3).^2 - Rhyd(r,6,0).^2); %at 10 eV
integrand369 = @(r) r.^2 .* Y(10,r) .* (Rhyd(r,6,1).^2 - Rhyd(r,6,0).^2); %at 10 eV

particleshift935 = integral(integrand935,0,50*RBOHR) %~1% error at m = 0, check vs mathematica...
particleshift369 = integral(integrand369,0,50*RBOHR) %~0.1% error at m = 0

masses = logspace(1,7,100)'; %10 eV-1e7 eV


ps369 = particleshift(r,dr,r.*Rhyd(r,6,1),r.*Rhyd(r,6,0),masses);
ps935 = particleshift(r,dr,r.*Rhyd(r,4,3),r.*Rhyd(r,6,0),masses);
ps369test = particleshift2(r,dr,(r.*Rhyd(r,6,1)).^2-(r.*Rhyd(r,6,0)).^2,masses);
ps935test = particleshift2(r,dr,(r.*Rhyd(r,4,3)).^2-(r.*Rhyd(r,6,0)).^2,masses);

%%% Compute Radial Density
%[r,dw,rR_6,~] = joon_tf_schrodinger(6,2,0,2000,5000);

%loglog(masses,abs(ps935),masses,abs(ps369))
%loglog(masses,1./(abs(ps935test - ps369test)))

% shifts935 = zeros(size(masses));
% shifts369 = zeros(size(masses));
% for k=1:numel(masses)
%     display(masses(k))
%     integrand935 = @(r) r.^2 .* Y(masses(k),r) .* (Rhyd(r,4,3).^2 - Rhyd(r,6,0).^2); %at zero mass for now
%     integrand369 = @(r) r.^2 .* Y(masses(k),r) .* (Rhyd(r,6,1).^2 - Rhyd(r,6,0).^2); %at zero mass for now
%     ps935(k) = integral(integrand935,0,100*RBOHR);
%     ps369(k) = integral(integrand369,0,1000*RBOHR);
% end
