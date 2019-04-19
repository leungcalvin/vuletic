isotopes = [168,170,172,174,176];
masses = [167.933897,169.9347618,171.9363815,173.9388621,175.9425717];
refidx = 3; %172;
mu = 1./masses - 1./masses(refidx);
m = 2.5;
b = 1e4;
g = 0;
dA = isotopes - isotopes(refidx);
noise = 1e-10; %in MHz
freqx = [-3,-2,0,2,3]*10^3; %MHz
shiftx= (freqx) ./ mu + randn(size(freqx)) * noise ./ mu;
shifty= shiftx * m + b + randn(size(shiftx)) .* noise;
nl = g .* dA .* masses .* masses(refidx) ./ (masses(refidx) - masses) ;
shifty = shifty + nl; %add king nonlinearity
freqy=shifty .* mu;freqy(refidx)=0;
freqx_err =noise * ones(size(freqx));
freqy_err =noise * ones(size(freqy));

[mw,bw,mwerr,bwerr,epsw,epswerr,rcs,pval] = kingtest(isotopes,masses,freqx,freqx_err,freqy,freqy_err);