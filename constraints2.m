RBOHR = 1;
Rhyd = @(r,n,l,netZ) sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2*netZ/(n.*RBOHR))^3)).*((2*netZ.*r)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*netZ.*r)./(n *RBOHR)) .*exp(-netZ.*r./(n*RBOHR));
Z = 70;
L = 100;
N = 2000;
[x,analytic70] = Rhyd2(2,4,Z,0,L,N);
[x,numeric70 ] = Rhyd2(2,10,Z,1,L,N);
dx = x(2) - x(1);

I0 = integral(@(r) Rhyd(r,4,2,2).^2 .* r.^2,0,Inf) %check that this integrates correctly
I1 = sum( analytic70(:,4).^2) * dx
I2 = sum( numeric70(:,4).^2 ) * dx

ax1=subplot(2,1,1);plot(x,analytic70.^2 ,x,Rhyd(x,4,2,70).^2 .* x.^2 );
ax2=subplot(2,1,2);plot(x,numeric70.^2,x,Rhyd(x,4,2,70).^2 .* x.^2  );
linkaxes([ax1,ax2],'xy')

n = 6;l = 0;
[x,rR_369g_test] = Rhyd2(l,n,70,0,10,10000);
R_369g = Rhyd(x,n,l,70);
figure;hold on;
ax1 =subplot(2,1,1);plot(x,rR_369g_test(:,n-l).^2,'green');
ax2 =subplot(2,1,2);plot(x,R_369g.^2 .* x.^2,'black');
linkaxes([ax1,ax2],'xy')
title(sprintf('n=%i,l=%i',n,l));

n = 6;l = 1;
[x,rR_369x_test] = Rhyd2(l,n,70,0,10,10000);
R_369x = Rhyd(x,n,l,70);
figure;hold on;
ax1 =subplot(2,1,1);plot(x,rR_369x_test(:,n-l).^2,'green');
ax2 =subplot(2,1,2);plot(x,R_369x.^2 .* x.^2,'black');
linkaxes([ax1,ax2],'xy')
title(sprintf('n=%i,l=%i',n,l));

n = 5;l = 2;
[x,rR_411x_test] = Rhyd2(l,n,70,0,10,10000);
R_411x = Rhyd(x,n,l,70);
figure;hold on;
ax3 =subplot(2,1,1);plot(x,rR_411x_test(:,n-l).^2,'green');
ax4 =subplot(2,1,2);plot(x,R_411x.^2 .* x.^2,'black');
linkaxes([ax3,ax4],'xy')
title(sprintf('n=%i,l=%i',n,l));

kingplotslope = 3.66;
