function fs = particleshift2(r,dr,drRsquared, a,c)
% rho(r) = rho0 * (1 + exp((r-c)/a))
%   c = 1.179157016706D-04 Bohr radii,
%   a = 9.890591365424D-06 Bohr radii;
%   rho0 doesn't really matter because we want ratios.
    radialDensity = (1 + exp((r - c)/a)).^-1;
    semilogx(r,radialDensity);
    fs = radialDensity .* drRsquared * dr';
end