function [x,rR] = Rhyd2(l_oam,nmodes,netZ,isNumerical,L,N)
    % l_oam = 0,1,2,3,...
    % nmodes =1,2,3,...
    % netZ = 2 for Yt+; 70 for Yt+ unscreened
    % isNumerical = 0 (analytical solution); = 1 (numerical solution)
    % L = box size in bohr radii. (default 100)
    % N = # of mesh points (default 10000). Need more for outer electron
    % orbitals.
    % L = 100;                   % Interval Length
    % N = 10000;                  % No of points
    x = linspace(0,L,N)';      % Coordinate vector
    dx = x(2) - x(1);          % Coordinate step
    rR = zeros(N,nmodes);
    if isNumerical == 0;
        'analytical'
        RBOHR = 1;
        %Rhyd = @(r,n,l,netZ) sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2*netZ/(n.*RBOHR))^3)).*((2*netZ.*r)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*netZ.*r)./(n *RBOHR)) .*exp(-netZ.*r./(n*RBOHR));
        for idx=(l_oam+1):1:(l_oam + nmodes)
            n = idx;
            l = l_oam;
            rR(:,idx-l_oam)= x.*sqrt((factorial(n - l - 1)./(2.*n.*factorial(n + l)).*(2*netZ/(n.*RBOHR))^3)).*((2*netZ.*x)./(n*RBOHR)).^l .* laguerreL(n-l-1,2.*l+1,(2*netZ.*x)./(n *RBOHR)) .*exp(-netZ.*x./(n*RBOHR));
            disp([n,l_oam,netZ]);
        end
        'done'
    else
        'numerical'

        % POTENTIAL, choose one or make your own
        %U = 1/2*100*x.^(2);    % quadratic harmonic oscillator potential
        %U = 1/2*x.^(4);       % quartic potential

        U = -netZ./x + l_oam*(l_oam+1) ./ x.^2; U(1) = 0; % Hydrogen, but with singularity at origin removed. Makes sense to zeroth order

        % Three-point finite-difference representation of Laplacian
        % using sparse matrices, where you save memory by only
        % storing non-zero matrix elements
        e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
        % Total Hamiltonian
        hbar = 1; m = 1;      % constants for Hamiltonian
        H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
        % Find lowest nmodes eigenvectors and eigenvalues of sparsematrix
        options.disp = 0;options.maxit = 20000;
        [rR,E] = eigs(H,nmodes,'sa',options);   % find eigs
        [E,ind] = sort(diag(E));% convert E to vector and sort low to high
        rR = rR(:,ind);           % rearrange corresponding eigenvectors
        rR = rR ./ sqrt(dx); % to match the analytical normalization
        'done'
    end

    end
