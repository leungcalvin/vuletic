L = 100;                      % Interval Length
N = 10000;                   % No of points
x = linspace(0,L,N)';      % Coordinate vector
dx = x(2) - x(1);           % Coordinate step
% POTENTIAL, choose one or make your own
%U = 1/2*100*x.^(2);    % quadratic harmonic oscillator potential
%U = 1/2*x.^(4);       % quartic potential
l_oam = 0;netZ = 2;nmodes = 10; n_principal = (1:1:nmodes)' + l_oam;
U = -netZ./x + l_oam*(l_oam+1) ./ x.^2; U(1) = 0;

% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
% Total Hamiltonian
hbar = 1; m = 1;      % constants for Hamiltonian
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
% Find lowest nmodes eigenvectors and eigenvalues of sparsematrix
options.disp = 0;options.maxit = 20000;options.tol = 1e-10;
[V,E] = eigs(H,nmodes,'sa',options);   % find eigs
[E,ind] = sort(diag(E));% convert E to vector and sort low to high
V = V(:,ind);           % rearrange corresponding eigenvectors
% Generate plot of lowest energy eigenvectors V(x) and U(x)
Usc = U*max(abs(V(:)))/max(abs(U));       % rescale U for plotting
figure;ax1=subplot(2,1,1);plot(x,V.^2);title(['Numerical (Hydrogenlike-ish) eigenfunctions, l=' string(l_oam)])
%plot(x,V.^2,x,Usc,'--k');          % plot V(x) and rescaled U(x)

% Add legend showing Energy of plotted V(x)
lgnd_str = [repmat('E = ',nmodes,1),num2str(E)];
legend(lgnd_str)                % place lengend string on plot
shg

RBOHR = 1;
Rhyd = @(r,n,l,netZ) sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2*netZ/(n.*RBOHR))^3)).*((2*netZ.*r)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*netZ.*r)./(n *RBOHR)) .*exp(-netZ.*r./(n*RBOHR));
ax2=subplot(2,1,2);title(['Analytical (Hydrogenlike) non-relativistic wavefunctions, l=' string(l_oam)]);hold on;
for idx=(l_oam+1):1:(l_oam + nmodes)
    plot(x,x.^2.*Rhyd(x,idx,l_oam,netZ).^2.*dx);
    disp([idx,l_oam,netZ]);
end
linkaxes([ax1,ax2],'xy')
