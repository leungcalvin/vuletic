function [r,dw,U,R] = joon_tf_schrodinger(l,Nlevel,rbi,rbf,Nr_ac)
% l = 2; % azimuthal angular momentum number
% Nlevel = 4; % number of levels
% rbi = 0; % start position of bin in a_0
% rbf = 2000; % end position of bin in a_0
rbcut = 5e-1; % cutpoint; below here grid is uniform with minumum binsize outside rcut
% Nr_ac = 5000; % number of r bins (and binpoints) after cutpoint rcut

% constants
alpha = 1/137.035999139; % fine structure constant
m_e = 1/alpha;

%% Define potential
% Thomas-Fermi potential for Yb+ (Z=70, n=2)
load('TF_Z=70_n=2.mat')
% Vfun = @(r) -Z*alpha./r;
Vfun = @(r) tfsol.V_TF(r);
% analytic solution
Efun = @(Z,n) -Z^2*alpha/2./n.^2;
Rfun = @(Z,n,l,r) sqrt(Z^3*(2/n^2)^2*factorial(n-l-1)/factorial(n+l))*exp(-Z*r/n).*(2*Z*r/n).^l.*laguerreL(n-l-1,2*l+1,2*Z*r/n);
Ufun = @(Z,n,l,r) r.*Rfun(Z,n,l,r);

%% calcuation parameters

if rbcut < rbi
    warning('rcut < ri.. rcut changed to ri')
    rbcut = rbi;
end
dlogw = log(rbf/rbcut)/Nr_ac;
dwmin = rbcut*(exp(dlogw)-1);
Nr_bc = round((rbcut - rbi)/dwmin); % binsize below cutpoint
dw_bc = (rbcut - rbi)/Nr_bc;
% r_bc = linspace(bi,bcut-dw_bc,Nr_bc);
% r_bc = linspace(bi+dw_bc/2,bcut-dw_bc/2,Nr_bc);
r_bc = linspace(rbi+dw_bc,rbcut,Nr_bc);
b_ac = exp(linspace(log(rbcut),log(rbf),Nr_ac+1));
dw_ac = diff(b_ac);
% r_ac = b_ac(1:end-1);
% r_ac = (b_ac(1:end-1)+b_ac(2:end))/2;
r_ac = b_ac(2:end);
Nr = Nr_bc + Nr_ac;
r = [r_bc,r_ac];

fprintf(' -- Nr = %u\n',Nr)

dr = diff(r); % distance between binpoints (dr(i) = r(i+1) - r(i))
dw = [dw_bc*ones(1,Nr_bc),dw_ac]; % width of each bin; useful for integration

%% solving equation
% second-order derivative -- non-uniform grid
w_l = +2./dr(1:end-1)./(dr(1:end-1)+dr(2:end)); % lower off-diagonal vector
w_d = -2./dr(1:end-1)./dr(2:end); % diagonal vector
w_u = +2./dr(2:end)./(dr(1:end-1)+dr(2:end)); % upper off-diagonal vector
% at boundary
w_l = [w_l,1/dr(end)^2,0];
w_d = [-2/dr(2)^2,w_d,-2/dr(end)^2];
w_u = [0,1/dr(1)^2,w_u];

DDM = spdiags([w_l.',w_d.',w_u.'],-1:1,Nr,Nr); % second-derivative matrix

% potential
V = Vfun(r);
Vmin = min(V);
Vin = V - Vmin; % to make potential always positive (negative potential results in unstability of solver somehow.)
VM = spdiags(Vin.',0,Nr,Nr); % potential matrix

% 1/r^2 term
IR2 = 1./r.^2;
IR2M = spdiags(IR2.',0,Nr,Nr); % 1/r^2 matrix

% matrix to be put in eigensolver
% radial schrodinger equation
S = -DDM + l*(l+1)*IR2M + 2*m_e*VM;

% eigensolver
eigoptions = struct;
eigoptions.issym = 0; % not symmetric for non-uniform grid
eigoptions.p = 30*Nlevel + 1;
eigoptions.isreal = 1;
opts.maxit = 500000; % default: 300
eigoptions.v0 = Ufun(2,l+1,l,r).'; % ground-state of hydrogen-like atom (Z=2)
% eigoptions.disp = 2;
tic
disp('---- eigensolver started -----')
[ev,ee,flag] = eigs(S,Nlevel,'sm',eigoptions);
fprintf('\t')
toc
disp('---- eigensolver ended -----')

if flag; warning('some eigenvalues are not converged.'); end

E = diag(ee)/2/m_e + Vmin; % energy levels
U = ev; % r*radial wavefunctions; U(:,n) is n-th wavefunction
U_analytic = zeros(Nr,Nlevel);
Z = 70;
for inx = 1:Nlevel
    norm2 = sum(dw.'.*U(:,inx).^2);
    U(:,inx) = U(:,inx)/sqrt(norm2); % normalization
    %U_analytic(:,inx) = r.* Rfun(Z,inx+l,l,r);
end
R = nan(Nr,Nlevel);
for inx = 1:Nlevel
    R(:,inx) = U(:,inx)./r.';
end

% order in energy
[E,I] = sort(E);
U = U(:,I);
R = R(:,I);
%U_analytic = U_analytic(:,I);

E
%% plot results
% match sign of each wavefunction
% [temp,rinx] = min(abs(r-0.1));% check near Bohr radius*0.1
% rinx = 1;
% r0 = r(rinx);
% for inx = 1:Nlevel
%     if U(rinx,inx)*Ufun(70,inx+l,l,r0) < 0 % Z=70
%         U(:,inx) = -U(:,inx);
%         R(:,inx) = -R(:,inx);
%     end
% end
% 
% cm = distinguishable_colors(Nlevel);
% wavefun_fig = figure;
% wavefun_fig.Position = [400,300,800,600];
% yyaxis left
% hold on
% 
% % n_draw = 6;
% n_draw = (l+1):(l+Nlevel);
% Inx_draw = [n_draw - l];
% 
% %calculation results
% h_cal = plot(r,R./r.'.^l,'-');
% for inx = 1:Nlevel
% %     h_cal(inx).Color = plotColor{inx}*0.6;
% %       plotColor{inx} = h_cal(inx).Color;
%       h_cal(inx).Color = cm(inx,:);
% %       h_cal(inx).Color = [1,0,0];
%       if ~ismember(inx,Inx_draw)
%           h_cal(inx).Visible = 'off';
%       end
% end
% % analytic solution
% Z = 70;
% for inx = 1:Nlevel
% %     if ismember(inx,inx_draw)
%     n = inx + l;
%     h_fun = plot(r,Rfun(Z,n,l,r)./r.^l);
%     h_fun.Color = [1,1,1]*0.5;
%     %     plotColor{inx} = h_fun.Color;
%     h_fun.LineStyle = '--';
%     h_fun.Marker = 'none';
%     if ~ismember(inx,Inx_draw)
%         h_fun.Visible = 'off';
%     end
% end
% Z = 2;
% for inx = 1:Nlevel
% %     if ismember(inx,inx_draw)
%     n = inx + l;
%     h_fun = plot(r,Ufun(Z,n,l,r)./r.^l);
%     h_fun.Color = [1,1,1]*0.5;
%     %     plotColor{inx} = h_fun.Color;
%     h_fun.LineStyle = ':';
%     h_fun.Marker = 'none';
%     if ~ismember(inx,Inx_draw)
%         h_fun.Visible = 'off';
%     end
% end
% ylabel('R(r)/r^l')
% 
% % potential
% yyaxis right
% h_V = plot(r,V);
% ylabel('\alphaV(r)/E_h')
% 
% yyaxis left
% 
% ax = wavefun_fig.CurrentAxes;
% ax.XScale = 'log';
% % ax.YScale = 'log';
% % xlim([bi,bf])
% xlim([rbi,1e2])
% xlabel('r/a_0')
% 
% title(sprintf(['Radial wavefunction of Yb+ w/ TF shielding (Z=%u, l=%u)\n',...
%     '[rbi,rbf]=[%g,%g], rbcut=%g, Nr\\_ac=%u'],Z,l,rbi,rbf,rbcut,Nr_ac));
% 
% Z = 70;
% n = (1:Nlevel).'+l;
% legendstr = compose('E_{%u%u} = %g (Z_{eff}=%g)',n,l,E,sqrt(-2*n.^2.*E/alpha));
% % legend(h_cal,legendstr{:},'position','best')
% legend(h_cal(Inx_draw),legendstr{Inx_draw});
%,'position','best')
%% analytic solution


% [Efun(2,(1:Nlevel)+l).',E,Efun(70,(1:Nlevel)+l).']/alpha/Z^2

end















