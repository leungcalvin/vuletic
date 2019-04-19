%%% potential considering helm dist of nuclear charge and Thomas-Fermi
%%% model for core electron potential

function V = V_total(r) % r in pm

%% potential due to Helm distribution of nuclear charge
Z = 70; % atomic number
A = 173; % mass number

alpha = 1/137.035999139;

a = 0.52; % fm
c = 1.23*A^(1/3) - 0.60; % fm
s = 0.9; % fm

r_N = sqrt(c^2 + 7/3*pi^2*a^2 - 5*s^2); % fm

%r: fm
V_H = @(r) -Z*alpha./(4*pi*r_N^3*r).*(sqrt(2*pi)*s.*((r.^2+r_N*r-2*r_N^2+2*s^2).*exp(-(r+r_N).^2/2/s^2)...
    - (r.^2-r_N*r-2*r_N^2+2*s^2).*exp(-(r-r_N).^2/2/s^2))...
    - pi*((r.^3-3*(r_N+s)*(r_N-s)*r-2*r_N^3).*erf((r+r_N)/sqrt(2)/s)...
    - (r.^3-3*(r_N+s)*(r_N-s)*r+2*r_N^3).*erf((r-r_N)/sqrt(2)/s)));

%% Thomas-Fermi model
Z = 70; % atomic number
n = 2; % ionization number (1 for singly ionized atom)
fname = sprintf('TF_Z=%u_n=%u.mat',Z,n);
if exist(fname,'file')
    load(fname);
else
    tfsol = TFsolver(Z,n);
    save(fname,'tfsol');
end

%% total potential -- pm^(-1); reduced compton wavelength of energy (E/hbar/c)
% r = logspace(-3,3,1000);

r0 = tfsol.r0;
r2x = tfsol.r2x;

V = zeros(1,length(r));
V(r < r0) = V_H(r(r<r0)*1e3)*1e3.*tfsol.chifun(r(r<r0)*r2x) - tfsol.n*alpha/r0;
V(r >= r0) = -tfsol.n*alpha./r(r>=r0);

% figure;
% loglog(r,-tfsol.n*alpha./r)
% hold on
% loglog(r,tfsol.V_TF(r))
% 
% loglog(r,-Z*alpha./r)
% loglog(r,V_H(r*1e3)*1e3)
% 
% loglog(r,V)
% 
% xlabel('r (pm)')
% ylabel('V/(hbar*c) (pm^{-1})')
% 
% ylim([-1e3,-5e-5])
% 
% legend('-n*\alpha/r',...
%     'Thomas-Fermi',...
%     '-Z*\alpha/r',...
%     'Helm',...
%     'Total',...
%     'location','best')

end