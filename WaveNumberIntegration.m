classdef WaveNumberIntegration < handle
    
    % Decide q's and k's and h
    % Calculate x's
    % Calculate p_i @ h
    % Calculate Phi
    % Calculate p_r @ z = h_1
    
    properties
        
        % Radius of transducer
        a
        % Distance from transducer to target
        z
        % Frequencies
        f
        % MultiLayerModel object
        model
        
    end
    
    properties (SetAccess = private)
        
        xi
        xo
        q
        theta
        
    end
    
    properties (Dependent = true)
        
        % Wavenumbers in front fluid
        kT
        % Wavenumbers in back fluid
        kR
        
    end
    
    methods
        
        function this = WaveNumberIntegration(xo, z, q, a, model)
            
            % Radius of transducer
            this.a = a;
            % Distance from source to target
            this.z = z;
            % Sine of angles and angles (rad)
            this.q = q;
            this.theta = asin(q);
            % Observation points (on the line with distance z from the
            % target
            this.xo = xo;
            % Material classes
            this.model = model;
            % x-positions where the incoming pressure field is calculated
            this.xi = z*tan(this.theta);
            
        end
        
        function [Pr, Pt] = doAll(this, f)
            % Combines the incoming pressure on the plate with the
            % reflection coefficients and propagates the pressure to the
            % observation points in xo.
            % 
            % [Pr, Pt] = doAll(this, f)
            % 
            % Input:
            % f: Frequencies
            % 
            % Output:
            % Pr: Reflected pressure (Number of observation points x Number
            % of frequencies)
            % Pt: Transmitted pressure (Number of observation points x Number
            % of frequencies)
            Pr = zeros(length(this.xo), length(f));
            Pt = zeros(length(this.xo), length(f));
            for i = 1:length(f);
                this.f = f(i);
                Phi = this.calcIncomingPressure(this.a);
                [V, W] = this.reflectionCoefficients();
                [Pr(:, i), Pt(:, i)] = this.reflectedPressure(Phi.', W, V);
            end
            
        end
        
        function Phi = calcIncomingPressure(this, a)
            % Calculates the incoming pressure from a plane piston of radius a.
            c = this.model.fluid(1).v;
            rho = this.model.fluid(1).density;
            U = 1/rho/c/2;
            f = this.f;
            q = this.q;
            theta = asin(q);
            
            z = this.z;
            xi = this.xi;
            % Calculate pressure p_i(x,z=h)
            Phi = zeros(length(theta), length(f));
            for i = 1:length(f)
                Phi(:, i) = planePistonPressure(xi, z, a, f(i), c, rho, U);
            end
            
            idz = f==0;
            if nnz(idz) > 0
                Phi(:, idz) = zeros(length(theta), 1);
            end
            
        end
        
        function [V, W] = reflectionCoefficients(this)
            % Calculates the reflection and transmission coefficients for a
            % fulid-solid-fluid system.
            %
            % [V, W] = reflectionCoefficients(this)
            %
            % Output:
            % V: Reflection coefficients
            % W: Transmission coefficients
            [V, W] = fluidSolidFluidReflectionCoefficient(this.f, this.theta, this.model);
        end
        
        function [Pr, Pt] = reflectedPressure(this, Phi, W, V)
            % Calculates the reflected and transmitted pressure at the
            % observation points xo at distance z from the target.
            % 
            % [Pr, Pt] = reflectedPressure(this, Phi, W, V)
            % 
            % Output:
            % Pr: Reflected pressure (Number of observation points x Number
            % of frequencies)
            % Pt: Reflected pressure (Number of observation points x Number
            % of frequencies)
            %
            % Input:
            % Phi: Incoming pressure
            % W: Transmission coefficient
            % V: Reflection coeficcient
            
            z = this.z;
            xo = this.xo;
            cR = this.model.fluid(1).v;
            cT = this.model.fluid(2).v;
            q = this.q;
            kR = this.kR;
            kT = this.kT;
            
            Pr = zeros(length(xo), length(kR));
            Pt = zeros(length(xo), length(kT));
            
            % Assume symmetric Phi, W and V and
            % that q = 0 is at the first element.
            W = [W(end:-1:2, :); W];
            V = [V(end:-1:2, :); V];
            q = [q(end:-1:2, :); q];
            Phi = [Phi(end:-1:2, :); Phi];
            for l = 1:length(xo)
                ER = exp((q.*xo(l) + sqrt(1-q.^2).*z)*1i*kR);
                ET = exp((q.*xo(l) + sqrt(1-q.^2).*z)*1i*kT);
                Iv = Phi.*V.*ER;
                Iw = Phi.*W.*ET;
                Pr(l, :) = trapz(q, Iv);
                Pt(l, :) = trapz(q, Iw);
            end
            
        end
        
        function k = get.kT(this)
            % Wavenumber in the front fluid
            c = this.model.fluid(1).v;
            k = 2*pi*this.f./c;
        end
        
        function k = get.kR(this)
            % Wavenumber in the back fluid
            c = this.model.fluid(2).v;
            k = 2*pi*this.f./c;
        end
        
    end
    
end