% Natural Units (+Gauss) radial Schrodinger solver
RBOHR=  5.29177210*10^(-11); %m
alpha = 1/137;
MASS = (alpha * RBOHR)^-1;
eps = RBOHR / 100;
omega = 1;
l = 0;
r = (0:eps:15*RBOHR)';
N = numel(r); %# of grid pts
s = 3; %# of energy levels to calculate
Lap = (-2*diag(ones(N,1),0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))./(eps.^2);
Lap(1,:) = 0; Lap(:,1) = 0; Lap(:,N) = 0; Lap(N,:) = 0;

netZ = 2;
V_eff = v_eff_unscreened2(r,l,netZ);V_eff(1) = 0;V_eff(N) = 0;

V_sho = 0.5 * MASS * omega.^2 * (r - mean(r)).^2;
% Schrodinger Eqn.
H = -Lap ./(2.* MASS) + spdiags(V_sho,0,N,N);
[u,lambda] = eigs(H,s);
%energies = diag(- lambda .*HBAR.^2./ (2.*MASS),0);
%R = repmat(1./r,[1,s]) .* u;
%Rhyd = @(r,n,l,netZ) sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2*netZ/(n.*RBOHR))^3)).*((2*netZ.*r)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*netZ.*r)./(n *RBOHR)) .*exp(-netZ.*r./(n*RBOHR));
for nr=1:s
  subplot(s,1,nr);
  plot(r/RBOHR,u(:,nr));hold on;
  %plot(r/RBOHR,r.^2.*Rhyd(r,l+nr+1,l).^2.*eps,'black');
end


