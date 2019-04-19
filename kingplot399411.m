% load shiftx,shifty, errx, erry, and fx = delta A / mu elementwise 
% so we can handle any pair of isotopes.

% 170 -> 172: -2044854.65843 +- 0.16555 kHz
% 172 -> 174: -1583068.32472 +- 0.24165 kHz
% 174 -> 176: -1509055.66623 +- 0.15257 kHz
% 170 -> 174: -3627922.42256 +- 0.15227 kHz
masses = zeros(1,176);
masses(168) = 167.933897;
masses(170) = 169.9347618;
masses(171) = 170.936583258;
masses(172) = 171.9363815;
masses(174) = 173.9388621;
masses(176) = 175.9425717;

A1 = [170,172,174,170];
A2 = [172,174,176,174];
freq411_MHz = [-2044854.65843,-1583068.32472,-1509055.66623,-3627922.42256].*1e-3;
err411_MHz = [0.16555,0.24165,0.15257,0.15227].*1e-3;

mu = (1./masses(A2) - 1./masses(A1)).*masses(172).^2;

%isotopes = [168,170,171,172,174,176];
absfreq399_MHz  = [1888.8, 1190.6, 941.3533, 531.11, 0.00, -508.89];
abserr399_MHz =   [0.11,   0.49,   0.250,    0.09,   0.09, 0.09];
freq399_MHz = [531.11-1190.6,0.00-531.11,-508.89-0.00,0.00-1190.6];
err399_MHz =  sqrt([0.09.^2 + 0.49.^2, 0.09.^2 + 0.250.^2, 0.09.^2 + 0.09.^2, 0.09.^2 + 0.49.^2]);

shiftx = freq411_MHz ./ mu;
shifty = freq399_MHz ./ mu;
errx = err411_MHz ./ mu;
erry = err399_MHz ./ mu;

fx = (A2 - A1) ./ mu;

[m_king,b_king,m_king_err,b_king_err,eps_king,eps_king_err, rcs_king,pval_king] = kingtest(shiftx,shifty,errx,erry,fx);
[m_lin,b_lin,m_lin_err,b_lin_err,  rcs_lin,pval_lin] = linearitytest(shiftx,shifty,errx,erry);
resid_lin = shifty - (b_lin + m_lin .* shiftx);
resid_lin_err = sqrt(b_lin_err.^2 + (m_lin_err .* shiftx).^2 + (m_lin.*errx).^2 + erry.^2);

resid_king = shifty - (b_king + m_king .* shiftx + eps_king .* fx);
resid_king_err = sqrt(b_king_err.^2 + (m_king_err .* shiftx).^2 + (m_king.*errx).^2 + erry.^2 +(eps_king_err.*fx).^2);

figure;
%subplot(3,1,3);hold on;
%stem(shiftx,resid_lin,'filled');
%ylabel('resid after king fit [freq-amu]');
%errorbar(shiftx,resid_lin,resid_lin_err,resid_lin_err,'o');

subplot(2,1,2);hold on;ylabel('resid [freq-amu]');
errorbar(shiftx,resid_king,resid_king_err,resid_king_err,'DisplayName','after king fit');
errorbar(shiftx,resid_lin,resid_lin_err,resid_lin_err,'DisplayName','after linear fit');
legend('-DynamicLegend');

subplot(2,1,1);hold on;nxsigma = 1;nysigma = 1;
scatter(shiftx,shifty);xlabel(sprintf('411 nm [MHz / amu]; err x %i',nxsigma));ylabel(sprintf('399 nm [MHz / amu]; err x %i',nysigma));
errorbar(shiftx,shifty,erry.*nysigma,erry.*nysigma,errx.*nxsigma,errx.*nxsigma,'o','DisplayName','Raw Data');
plot(shiftx,b_lin+m_lin.*shiftx,'DisplayName',sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e)\n[rcs =%0.2e,pval=%0.2e]',m_lin,m_lin_err,b_lin,b_lin_err,rcs_lin,pval_lin),'color','red');
plot(shiftx,b_king+m_king.*shiftx+eps_king.*fx,'DisplayName',sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e) + (%0.6e +/- %0.6e)f(x)\n[rcs =%0.2e,pval=%0.4e]',m_king,m_king_err,b_king,b_king_err,eps_king,eps_king_err,rcs_king,pval_king),'color','green');
%sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e)\n[rcs =%0.2e,pval=%0.2e]',mw,mwerr,bw,bwerr,rcs,pval);
legend('-DynamicLegend')
title('King Plot + Linear Fits');
