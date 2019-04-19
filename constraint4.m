% Looks at Li I levels as in Example 1 in the Grasp2K Manual
% example1 -> example4 for Yb nuclei with same electronic structure
nucleus = '(Z=70)'; N = 253 % grid points
%nucleus = '(Z=3)'; N = 333; % grid points

grasp1s = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example4/plot_1s.dat', 1, N));
r1s = grasp1s(:,1);large_1s = grasp1s(:,2);small_1s = grasp1s(:,3);density_1s = large_1s.^2 + small_1s.^2;

grasp2p12 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example4/plot_2p12.dat', 1, N)); %first Ngrid points are all the same...
r2p12 = grasp2p12(:,1);large_2p12 = grasp2p12(:,2);small_2p12 = grasp2p12(:,3);density_2p12 = large_2p12.^2 + small_2p12.^2;

grasp2p32 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example4/plot_2p32.dat', 1, N));
r2p32 = grasp2p32(:,1);large_2p32 = grasp2p32(:,2);small_2p32 = grasp2p32(:,3);density_2p32 = large_2p32.^2 + small_2p32.^2;


assert(sum((r1s - r2p12).^2) == 0);
assert(sum((r1s - r2p32).^2) == 0);
r = r1s; dw = [0; diff(r)]; %putting integration measure dw = 0 at the origin, be careful when looking at masses scales heavier than (r(2) - r(1)) bohr radii

loglog(r,large_1s,'DisplayName','1s');hold on;
plot(r,large_1s(2).*(r./r(2)),'DisplayName','r^1');

plot(r,large_2p32,'DisplayName','1s');
plot(r,large_2p32(2).*(r./r(2)).^2,'DisplayName','r^2');

% Transition A: 1s -> 2p1/2
% Transition B: 1s -> 2p3/2
drRsquaredA = 1 .* (density_2p12 - density_1s); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)
drRsquaredB = 1 .* (density_2p32 - density_1s); % 

n_m = 1000;
particlemasses = logspace(1,log10(3e7),n_m)'; %10 eV-1e7 eV
psA= particleshift2(r',dw',drRsquaredA',particlemasses);
psB= particleshift2(r',dw',drRsquaredB',particlemasses);

%% A/B comparison
kingslope = 1;power = 4;
constraint = 1./(abs(kingslope*psA - psB));
loglog(particlemasses,constraint,'DisplayName',['GRASP2K wavefunctions ', nucleus]);hold on;
loglog(particlemasses,particlemasses.^power .* constraint(numel(constraint))./particlemasses(n_m).^power,'DisplayName',['powerlaw exp='+string(power)])

title(['kingslope = ', string(kingslope)]);
%axis([min(particlemasses),max(particlemasses),1e-18,1e-6]);
legend('-DynamicLegend');

axis([min(particlemasses),max(particlemasses),1e-18,1e-0]);

