%% One Electron Overlap Integrals
N=336;%N=353;
grasp1S0 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5/plot_6s_1S0.dat', 1, N));
r1S0 = grasp1S0(:,1);large_1S0 = grasp1S0(:,2);small_1S0 = grasp1S0(:,3);density_1S0 = large_1S0.^2 + small_1S0.^2;

N=336;%N=352;
grasp1P0 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5/plot_6p_1P0.dat', 1, N));
r1P0 = grasp1P0(:,1);large_1P0 = grasp1P0(:,2);small_1P0 = grasp1P0(:,3);density_1P0 = large_1P0.^2 + small_1P0.^2;

N=336;
grasp1P1 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5/plot_6p_1P1.dat', 1, N));
r1P1 = grasp1P1(:,1);large_1P1 = grasp1P1(:,2);small_1P1 = grasp1P1(:,3);density_1P1 = large_1P1.^2 + small_1P1.^2;

rneutral = r1S0; dwneutral = [0;diff(r1S0)];
wfcheck = 0;
if wfcheck %check that large component wavefunction satisfies P(r) ~ r*R(r) ~ r^(l+1) for small r
    %powercheck(rneutral,large_1S0,2,1);
    %powercheck(rneutral,large_1P0,2,2);
end

drRsquared399 = 1 .* (density_1P1 - density_1S0); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)


%% 435/411 transition
N=417;%N=336
grasp2S12 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_6s_2S12.dat', 1, N));
r2S12 = grasp2S12(:,1);large_2S12 = grasp2S12(:,2);small_2S12 = grasp2S12(:,3);density_2S12 = large_2S12.^2 + small_2S12.^2;
%powercheck(r2S12,large_2S12,2,1);

N=417;%N=340;
grasp2D32 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_5d_2D32.dat', 1, N));
r2D32 = grasp2D32(:,1);large_2D32 = grasp2D32(:,2);small_2D32 = grasp2D32(:,3);density_2D32 = large_2D32.^2 + small_2D32.^2;
%powercheck(r2D32,large_2D32,2,3);

N=417;%N=355;
grasp2D52 = table2array(import_rwf('/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_5d_2D52.dat', 1, N));
r2D52 = grasp2D52(:,1);large_2D52 = grasp2D52(:,2);small_2D52 = grasp2D52(:,3);density_2D52 = large_2D52.^2 + small_2D52.^2;
%powercheck(r2D52,large_2D52,2,3);

rion = r2S12; dwion = [0;diff(r2S12)];

drRsquared411 = 1 .* (density_2D52 - density_2S12); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)
drRsquared435 = 1 .* (density_2D32 - density_2S12); % GRASP2K spits out the upper and lower components of u(r) = r*R(r)

n_m = 1000;
particlemasses = logspace(1,log10(3e9),n_m)'; %10 eV-1e7 eV
ps399= particleshift2(rneutral',dwneutral',drRsquared399',particlemasses);
ps435= particleshift2(rion',dwion',drRsquared435',particlemasses);
ps411= particleshift2(rion',dwion',drRsquared411',particlemasses);

%% A/B comparison
kingslope = 1;
nl_eps =4.51e6 ; %Hz;
curve435vs411=constraintcurve(ps411,ps435,1,1);
loglog(particlemasses,curve435vs411,'DisplayName',['GRASP2K wavefunctions (x=411,y=435,F_y/F_x=0.27)']);hold on;


%Always use 411 on x axis of king plot
for idx=1:40
    curve399vs411=constraintcurve(ps411,ps399,1,1);
    loglog(particlemasses,curve399vs411,'DisplayName',['(x=411,y=399,F_y/F_x='+string(0.272 + 0.05*idx)]);hold on;
end



