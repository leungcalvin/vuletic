%% Transition X: 397 nm 2S1/2->2P1/2
N=321;%N=333;
grasp2S12 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_4s_2S12.dat', 1, N));
r2S12 = grasp2S12(:,1);large_2S12 = grasp2S12(:,2);small_2S12 = grasp2S12(:,3);density_2S12 = large_2S12.^2 + small_2S12.^2;

N=321;%N=321;
grasp2P12 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_4p_2P12.dat', 1, N));
r2P12 = grasp2P12(:,1);large_2P12 = grasp2P12(:,2);small_2P12 = grasp2P12(:,3);density_2P12 = large_2P12.^2 + small_2P12.^2;

N=321;%N=321;
grasp2P32 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_4p_2P32.dat', 1, N));
r2P32 = grasp2P32(:,1);large_2P32 = grasp2P32(:,2);small_2P32 = grasp2P32(:,3);density_2P32 = large_2P32.^2 + small_2P32.^2;

rneutral = r2S12; dwneutral = [0;diff(r2S12)];
plot(rneutral,(large_2S12),'DisplayName','2S12');hold on;
%plot(r, large_2S12(2) .* (r ./ r(2)),'DisplayName','r^1'); %wavefunctions are r*R(r) ~ r^(l+1) for small r

plot(rneutral,(large_2P32),'Displayname','2P32');hold on;
%plot(r, large_2P32(2) .* (r ./ r(2)).^2,'DisplayName','r^2');

plot(rneutral,(large_2P12),'Displayname','2P12');hold on;
%plot(r, large_2P12(2) .* (r ./ r(2)).^2,'DisplayName','r^2');
legend('-DynamicLegend')

drRsquaredX = 1 .* (density_2P12 - density_2S12); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)
% plot(r,density_2P32,'DisplayName','2P32 density');hold on;
% plot(r,density_2S12,'DisplayName','2S12 density');

%% 866 nm 2P1/2->2D3/2
N=321;
grasp2P12 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_4p_2P12.dat', 1, N));
r2P12 = grasp2P12(:,1);large_2P12 = grasp2P12(:,2);small_2P12 = grasp2P12(:,3);density_2P12 = large_2P12.^2 + small_2P12.^2;

grasp2D32 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_3d_2D32.dat', 1, N));
r2D32 = grasp2D32(:,1);large_2D32 = grasp2D32(:,2);small_2D32 = grasp2D32(:,3);density_2D32 = large_2D32.^2 + small_2D32.^2;

grasp2D52 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/calcium/plot_3d_2D52.dat', 1, N));
r2D52 = grasp2D52(:,1);large_2D52 = grasp2D52(:,2);small_2D52 = grasp2D52(:,3);density_2D52 = large_2D52.^2 + small_2D52.^2;

rion = r2P12; dwion = [0;diff(r2P12)];
plot(rion,(large_2S12),'DisplayName','2S12');hold on;
%plot(r, large_2S12(2) .* (r ./ r(2)),'DisplayName','r^1'); %wavefunctions are r*R(r) ~ r^(l+1) for small r

plot(rion,(large_2D32),'Displayname','2D32');hold on;
%plot(r, large_2D32(2) .* (r ./ r(2)).^2,'DisplayName','r^2');

plot(rion,(large_2D52),'Displayname','2D32');hold on;
%plot(r, large_2D32(2) .* (r ./ r(2)).^2,'DisplayName','r^2');

legend('-DynamicLegend')

drRsquaredY = 1 .* (density_2D32 - density_2P12); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)

n_m = 1000;
particlemasses = logspace(1,log10(3e9),n_m)'; %10 eV-1e7 eV
psY= particleshift2(rion',dwion',drRsquaredY',particlemasses);
psX= particleshift2(rneutral',dwneutral',drRsquaredX',particlemasses);

%% A/B comparison
kingslope = -0.3112136266855926;%slope positive or negative?
power = 2;nl_eps = 0.1e6; %Hz
constraint = nl_eps./(abs(psY - kingslope .* psX));
loglog(particlemasses,constraint,'DisplayName',['GRASP2K wavefunctions 397/866']);hold on;
loglog(particlemasses,particlemasses.^power .* constraint(numel(constraint))./particlemasses(n_m).^power,'DisplayName',['powerlaw exp='+string(power)])

title(['kingslope = ', string(kingslope)]);
axis([min(particlemasses),max(particlemasses),1e-1*constraint(1),1e1*constraint(n_m)]);
legend('-DynamicLegend');

