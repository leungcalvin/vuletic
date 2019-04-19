function ps = particleshift2(r,dr,drRsquared, masses)
    % r,dr, rR_ground, rR_excited: (1 x n) row vectors, 
    % masses as (m x 1) column vector.
    % hbar = 1, c = 1
    % remember that F(r), G(r) are relativistic analogues of u(r) = r*R(r);
    N = numel(r);
    if r(1) == 0
        r = r(:,2:N);
        dr = dr(:,2:N);
        drRsquared = drRsquared(:,2:N);
    end
    SL = 3e8; %m/s
    HBAR = 1.054e-34; %J-s
    RBOHR=  5.29177210*10^(-11); %m
    EVperJ = 1.602*10^(-19); %1 eV / 1 Joule
%     size(masses)
%     size(drRsquared) %drRsquared = (n x 1)
%     size(r)
%     size(dr)
%     size(drRsquared)
    %disp(drRsquared)
    ps = (2.*pi).^-1.*SL.*RBOHR.^(-1).*exp(-masses *r ./ 3728.939) * (dr.*drRsquared./(4*pi.*r))';
%     figure;
%     plot(r,exp(-EVperJ.*masses *r ./ (HBAR* SL)));
%     figure;
%     plot(r,(dr.*drRsquared./(4*pi.*r))');
%     %ps = (2.*pi).^-1.*SL.*exp(-EVperJ.*masses *r ./ (HBAR* SL)) * (drRsquared./(4*pi.*r))' .* dr;% masses = (mx1), r = (1xn), Y = (m x n)
end