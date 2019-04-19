%% One Electron Overlap Integrals V2
%% Ground state has occupation number of -1, excited state has occupation number of +1.
%% One Electron Densities

[grid435,density435] = getGRASP2KDensity(417,{'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_6s_2S12.dat',
                                              '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_5d_2D32.dat'},[-1;1]);dw435 = midbin(grid435);

                                          
[grid411,density411] = getGRASP2KDensity(417,{'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_6s_2S12.dat',
                                              '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/plot_5d_2D52.dat'},[-1;1]);dw411 = midbin(grid411);
                                    
[grid399,density399] = getGRASP2KDensity(336,{'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_6s_1S0.dat';
                                              '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_6p_1P0.dat'},[-1;1]);dw399 = [0;diff(grid399)];

%% 399 Multi-electron density
filelist_399 = {'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_5p_2P12.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_5p_2P32.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_5s_1S0.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_6p_1P0.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_6p_1P1.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/ground/plot_6s_1S0.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_5p_2P12.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_5p_2P32.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_5s_1S0.dat';
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_6p_1P0.dat'; %399 excited state
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_6p_1P1.dat'; %triplet state
'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/example5test/excited/plot_6s_1S0.dat';
};
occupation_399= [-2;-4;-2;0;0;-2; 2;4;2;1;0;1];
[grid399_multi,density399_multi] = getGRASP2KDensity(332,filelist_399,occupation_399);dw399_multi = midbin(grid399_multi);
%% 411 Multi-electron density
filelist_411 = {'/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5s_2S12.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5p_2S12.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5p-_2S12.dat';    
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_6s_2S12.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5s_2D52.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5p-_2D52.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5p_2D52.dat';
    '/home/calvin/Documents/vuletic/grasp2K_v1_1/manual/ion/multielectron/411/rwfn_5d_2D52.dat';
    }
occupation_411 = [-2, -4,-2, -1, 2,2,4,1];

[grid411_multi,density411_multi] = getGRASP2KDensity(154,filelist_411,occupation_411);dw411_multi = midbin(grid411_multi);
%% Particle Shifts
n_m = 1000;
particlemasses = logspace(1,log10(3e7),n_m)'; %10 eV-1e7 eV
ps399= particleshift2(grid399',dw399',density399',particlemasses);
ps399_multi=particleshift2(grid399_multi',dw399_multi',density399_multi',particlemasses);
ps435= particleshift2(grid435',dw435',density435',particlemasses);
ps411= particleshift2(grid411',dw411',density411',particlemasses);
ps411_multi= particleshift2(grid411_multi',dw411_multi',density411_multi',particlemasses);

%% F Ratios
fs399 = fieldshift(grid399',dw399',density399',9.89059e-6,1.1791570e-4); % see 6s_6p_DF.sum for nuclear params
fs399_multi= fieldshift(grid399_multi',dw399_multi',density399_multi',9.89059e-6,1.1791570e-4); % see 6s_6p_DF.sum for nuclear params
fs411 = fieldshift(grid411',dw411',density411',9.89059e-6,1.1791570e-4); % see 6s_6p_DF.sum for nuclear params
fs411_multi= fieldshift(grid411_multi',dw411_multi',density411_multi',9.89059e-6,1.1791570e-4); % see 6s_6p_DF.sum for nuclear params
fs435 = fieldshift(grid435',dw435',density435',9.89059e-6,1.1791570e-4); % see 6s_6p_DF.sum for nuclear params
eps_king = 4.5145e6; % for 411 399 comparison
m_king = 0.273; % for 411 399 experimental value

loglog(particlemasses,constraintcurve(ps411,ps399,fs399/fs411,1),'DisplayName','(x=411,y=399)');hold on;
loglog(particlemasses,constraintcurve(ps435,ps399,fs399/fs435,1),'DisplayName','(x=435,y=399)');
loglog(particlemasses,constraintcurve(ps411,ps399_multi,fs399_multi/fs411,1),'DisplayName','(x=411,y=399m)');
%loglog(particlemasses,constraintcurve(ps411_multi,ps399_multi,fs399_multi/fs411_multi,1),'DisplayName','(x=411m,y=399m)');
loglog(particlemasses,constraintcurve(ps435,ps399_multi,fs399_multi/fs435,1),'DisplayName','(x=435,y=399m)');
loglog(particlemasses,constraintcurve(ps435,ps411,fs411/fs435,1),'DisplayName','(x=435,y=411)');
atomki; %plot atomki anomaly from 1602.04822
legend('-DynamicLegend');
title('eps = 1 Hz sensitivity');
axis([1e1,max(particlemasses),1e-16,1e-6]);
%the densities here are defined as r^2 R(r)^2.

