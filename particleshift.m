function ps = particleshift(r,dr,rR_ground,rR_excited, masses)
    % r,dr, rR_ground, rR_excited: (1 x n) row vectors, 
    % masses as (m x 1) column vector.
    N = numel(r);
    r = r(:,2:N);
    dr = dr(:,2:N);
    rR_ground = rR_ground(:,2:N); %omit the point at r = 0
    rR_excited= rR_excited(:,2:N);
    SL = 3e8; %m/s
    HBAR = 1.054e-34; %J-s
    %RBOHR=  5.29177210*10^(-11); %m
    EVperJ = 1.602*10^(-19); %1 eV / 1 Joule
    drRsquared = rR_ground.^2 - rR_excited.^2;
%     size(masses)
%     size(drRsquared) %drRsquared = (n x 1)
%     size(r)
%     size(dr)
%     size(drRsquared)
    ps = (2.*pi).^-1.*SL.*exp(-EVperJ.*masses *r ./ (HBAR* SL)) * (dr.*drRsquared./(4*pi.*r))';
    %ps = (2.*pi).^-1.*SL.*exp(-EVperJ.*masses *r ./ (HBAR* SL)) * (drRsquared./(4*pi.*r))' .* dr;% masses = (mx1), r = (1xn), Y = (m x n)
end