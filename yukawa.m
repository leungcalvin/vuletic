function Y = yukawa(r,m)
    assert(length(m) == 1);
    assert(length(l) == 1);
    SL = 3e8; %m/s
    HBAR = 1.054e-34; %J-s
    EVperJ = 1.602*10^(-19); %1 eV / 1 Joule
    Y = (2*pi)*SL*exp(-m*EVperJ*r ./ (HBAR* SL)) ./ (4*pi*r)
end