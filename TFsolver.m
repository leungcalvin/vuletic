%% TF2_bvp4c.m being class
classdef TFsolver < handle
    properties (SetAccess=immutable)
        Z % atomic number
        n % ionization number (1 for singly ionized atom)
    end
    properties (SetAccess=private)
        x0 % correction point
        r2x
        r0
    end
    properties
        chifun % chi as a function w.r.t x: normalized r with x0
%         V_TF % potential due to core electrons as a function w.r.t r in pm
    end
    properties(Constant)
        alpha = 1/137.035999139; % fine structure constant
        E_e = 0.5109989461; % MeV; Mass energy of electron
        hbarc = 0.19732697; % eV um = Mev pm; hbar*c
    end
    
    methods
        %%% Constructor
        function obj = TFsolver(Z,n)
            x0 = fzero(@(x0) TFsolver.bvpsolve_Dyx0_res(x0,Z,n),7);
            % error = bvpsolve_Dyx0_res(x0_sol);
            
            sol = TFsolver.bvpsolve(x0);
            chifun = @(x) TFsolver.deval_first(sol,x);
            
            obj.Z = Z;
            obj.n = n;
            obj.chifun = chifun;
            obj.x0 = x0;
            
            m_e = obj.E_e/obj.hbarc; % pm^(-1); inverse of reduced compton wavelength of electron (1/lambda-bar)
            
            r2x = 4*(2*obj.Z/9/pi^2)^(1/3)*m_e*obj.alpha;
            r0 = obj.x0/r2x; % pm
            
            obj.r2x = r2x;
            obj.r0 = r0;
        end
        
        % potential due to core electrons as a function w.r.t r in pm
        % unit of potential energy: pm^(-1); reduced compton wavelength of energy (E/hbar/c)
        function V = V_TF(obj,r)
            V = zeros(1,length(r));
            V(r < obj.r0) = -obj.Z*obj.alpha./r(r<obj.r0).*obj.chifun(r(r<obj.r0)*obj.r2x) - obj.n*obj.alpha/obj.r0;
            V(r >= obj.r0) = -obj.n*obj.alpha./r(r>=obj.r0);
        end
    end
    
    methods(Static)
        function y1 = deval_first(sol,x)
            y = real(deval(sol,x));
            y1 = y(1,:);
        end
        
        function Dyx0 = bvpsolve_Dyx0_res(x0,Z,n)
            [sol,Dyx0] = TFsolver.bvpsolve(x0);
            Dyx0 = real(Dyx0) + n/Z;
        end
        
        function [sol,Dyx0] = bvpsolve(x0)
            
            % x0 = 10;
            xrange = [0,x0];
            
            solinit = bvpinit(linspace(xrange(1),xrange(2),5),[1 0]);
            sol = bvp4c(@TFsolver.twoode,@TFsolver.twobc,solinit);
            
            % x = linspace(xrange(1),xrange(2),100);
            % y = deval(sol,x);
            
            y = deval(sol,x0);
            
            Dyx0 = y(2,end);
            
            % figure;
            % plot(x,y(1,:),'-')
        end
        
        function dydx = twoode(x,y)
            dydx = [ y(2); sqrt(y(1)^3/(x+1e-10)) ]; % x+1e-10 instead x to avoid singularity
        end
        
        function res = twobc(ya,yb)
            res = [ya(1) - 1; yb(1)];
        end
    end
end

