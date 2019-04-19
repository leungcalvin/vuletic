function [shiftx,shifty,errx,erry] = kingplot3(isotopes,masses, fx,dfx,fy,dfy)
% choose reference isotope
[~,refidx] = max(isotopes == 172)
% name each point on the plot
mu = (1./masses - 1/masses(refidx))*masses(refidx).^2;
dA = isotopes - isotopes(refidx);
names = string(isotopes) + '-' + string(isotopes(refidx));
shiftx = (fx - fx(refidx)) ./ mu;
shifty = (fy - fy(refidx)) ./ mu;
errx = sqrt(dfx.^2 + dfx(refidx).^2)./mu;
erry = sqrt(dfy.^2 + dfy(refidx).^2)./mu;
fx = dA .* masses ./( masses(refidx) .* (masses(refidx) - masses));

% remove nans & reference isotope...
keepers = ~boolean(isnan(shiftx + shifty)); %
shiftx = shiftx(keepers);
shifty = shifty(keepers);
errx = errx(keepers);
erry = erry(keepers);
names = names(keepers);

fx = fx(keepers);
%
[m_king,b_king,m_king_err,b_king_err,eps_king,eps_king_err, rcs_king,pval_king] = kingtest(shiftx,shifty,errx,erry,fx)
[m_lin,b_lin,m_lin_err,b_lin_err,  rcs_lin,pval_lin] = linearitytest(shiftx,shifty,errx,erry)
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
errorbar(shiftx,resid_king,resid_king_err,resid_king_err,'DisplayName','after king fit','color','green');
errorbar(shiftx,resid_lin,resid_lin_err,resid_lin_err,'DisplayName','after linear fit','color','red');
legend('-DynamicLegend');

subplot(2,1,1);hold on;
scatter(shiftx,shifty);
errorbar(shiftx,shifty,erry,erry,errx,errx,'o','DisplayName','Raw Data');
plot(shiftx,b_lin+m_lin.*shiftx,'DisplayName',sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e)\n[rcs =%0.2e,pval=%0.2e]',m_lin,m_lin_err,b_lin,b_lin_err,rcs_lin,pval_lin),'color','red');
plot(shiftx,b_king+m_king.*shiftx+eps_king.*fx,'DisplayName',sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e) + (%0.6e +/- %0.6e)f(x)\n[rcs =%0.2e,pval=%0.4e]',m_king,m_king_err,b_king,b_king_err,eps_king,eps_king_err,rcs_king,pval_king),'color','green');
%sprintf('y = (%0.6e +/- %0.6e)x + (%0.6e +/- %0.6e)\n[rcs =%0.2e,pval=%0.2e]',mw,mwerr,bw,bwerr,rcs,pval);
legend('-DynamicLegend')
title('King Plot + Linear Fits');
text(shiftx, shifty, names, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');


end

