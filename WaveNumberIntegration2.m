classdef WaveNumberIntegration2 < handle
    
    properties
        
        % Distance from transducer (TX) to target
        h
        
    end
    
    properties (Dependent = true)
        
        % Model
        model
        % Trasducer (TX) radius
        a
        % Reflection and transmission coefficient
        R
        T
        % Incoming pressure field (wavenumber domain)
        Phi
        % Horizontal and vertical wave numbers for all frequencies and
        % angles (matrix of size nq x nf, nq = length(q), nf = length(f))
        kz
        kx
        % Sine of angle of the plane wave
        q
        % Frequencies in the model
        f
        % Number of plane waves
        nq
        
    end

    properties (Access = protected)
        
        % Variables for keeping track of whether R, T or Phi needs to be
        % updated.
        updateRT
        updatePhi
        
        % Variables for storing calculated reflection and transmission
        % coefficient and the incoming pressure (k-domain)
        R_
        T_
        Phi_
        a_
        model_
        f_
        nq_
        
    end
    
    methods
   
        function this = WaveNumberIntegration2(a, h, f, model, nq)
            
            % Set transducer (TX) radius
            this.a = a;
            % Set distance from TX to target
            this.h = h;
            % Set frequencies
            this.f = f;
            this.model = model;
            % Set the number of sampling points in k-space
            this.nq = nq;
            
        end
        
        function P = calculateReflectedPressure(this, x, z, R)
            % Calculates the pressure in each pair of points {x(i), z(i)} in frequency domain.
            % z is the distance from the top of the plate.
            % 
            % Input:
            % x - x-positions
            % z - z-positions
            % R - (Optional) Either a reflection coefficient (default) or a
            %     transmission coefficient.
            if nargin == 3
                R = this.R;
            end
            % Frequency
            f = this.f();
                        
            rho = this.model.fluid(1).density;
            Kx = this.kx;
            Kz = this.kz;
            
            q = this.q;
            
            V = this.Phi;
            % Add the distance from the transducer (TX) to the z position
            % to take into account the propagation from TX to target.
            z = this.h + z;
            
            P = propagateReflectedWave(x, z, f, q, Kx, Kz, rho, V, R);
        end
        
        function P = calculateTransmittedPressure(this, x, z)
            % Calculates the pressure in each pair of points {x(i), z(i)} in frequency domain.
            % z is the distance from the bottom of the plate.
            T = this.T;
            P = this.calculateReflectedPressure(x, z, T);
        end
        
        function [x, t] = calculateReflectedPressureRx(this, a, d, nx, X, fs, fx)
            % Calculates the time signal experienced on the transducer
            % surface.
            %
            % x - Time signal
            % t - Times (s)
            xx = linspace(0, a, nx);
            zz = d;
            % Get the pressure (matrix nx x nf)
            P = this.calculateReflectedPressure(xx, zz);
            
            % Convolve with the pulse
            nfft = size(P, 2);
            c = zeros(nx, nfft);
            for i = 1:nx
                [c(i, :), t] = convTimeFreq(P(i, :), fs, nfft, X, fx, true);
            end
            
            % Integrate over the transducer surface. Assume circular
            % symmetric???
            x = trapz(xx, 2*pi*xx*c, 1);
        end
        
        function [x, t] = calculateTransmittedPressureRx(this, d, a)
            
        end
        
        function kx = get.kx(this)
            % Returns a matrix of horizontal wave numbers in the front
            % fluid
            w = 2*pi*this.f;
            q = this.q(:);
            nq = length(q);
            nf = length(w);
            
            k = w./this.model.fluid(1).v;
            K = repmat(k, nq, 1);
            Q = repmat(q, 1, nf);
            kx = K.*Q;
        end
        
        function kz = get.kz(this)
            % Returns a matrix of vertical wave numbers in the front
            % fluid
            w = 2*pi*this.f;
            cq = sqrt((1 - this.q(:).^2));
            nq = length(cq);
            nf = length(w);
            
            k = w./this.model.fluid(1).v;
            K = repmat(k, nq, 1);
            CQ = repmat(cq, 1, nf);
            kz = K.*CQ;
        end
        
        function Phi = get.Phi(this)
            % Updates the reflection and transmission coeffecient if updateRT
            % is true. Otherwise returns the variable T_.
            
            if this.updatePhi
                a = this.a;
                c = this.model.fluid(1).v;
                f = this.f;
                q = this.q;
                nf = length(f);
                Phi = zeros(this.nq, nf);
                Kx = this.kx;
                for i = 1:nf
                    kx = Kx(:, i);
                    Phi(:, i) = angularPlaneWaveSpectrumPiston(a, c, q, f);
                    this.Phi_ = Phi;
                    this.updatePhi = false;
                end
            else
                Phi = this.Phi_;
            end
            
        end
        
        function R = get.R(this)
            % Updates the reflection and transmission coeffecient if
            % updateRT is true. Otherwise returns the variable R_.
            if this.updateRT
                [R, T] = fluidSolidFluidReflectionCoefficient(this.f, asin(this.q), this.model);
                this.R_ = R;
                this.T_ = T;
                this.updateRT = false;
            else
                R = this.R_;
            end
        end
        
        function T = get.T(this)
            % Updates the reflection and transmission coeffecient if updateRT
            % is true. Otherwise returns the variable T_.
            if this.updateRT
                [R, T] = fluidSolidFluidReflectionCoefficient(this.f, asin(this.q), this.model);
                this.R_ = R;
                this.T_ = T;
                this.updateRT = false;
            else
                T = this.T_;
            end
        end
        
        function set.model(this, model)
            % Update model parameters
            
            % Maybe do some checking on model here.
            this.model_ = model;
            this.updateRT = true;
        end
        
        function model = get.model(this)
            % Returns the layered model
            model = this.model_;
        end
        
        function set.a(this, val)
            % Set the transducer (TX) radius (m)
            this.a_ = val;
            this.updatePhi = true;
        end
        
        function val = get.a(this)
            % Returns the transducer (TX) radius (m)
            val = this.a_;
        end
        
        function val = get.q(this)
            % Get's the sine of angle of the plane waves
            val = linspace(-0.99, 0.99, this.nq);
        end
        
        function set.f(this, val)
            % Sets the frequencies in the model and sets the updateRT and
            % updatePhi variables to true.
            this.f_ = val;
            this.updateRT = true;
            this.updatePhi = true;
        end
        
        function val = get.f(this)
            % Returns the frequencies in the model
            val = this.f_;
        end

        function set.nq(this, val)
            % Sets the frequencies in the model and sets the updateRT and
            % updatePhi variables to true.
            this.nq_ = val;
            this.updateRT = true;
            this.updatePhi = true;
        end
        
        function val = get.nq(this)
            % Returns the frequencies in the model
            val = this.nq_;
        end
        
    end
    
end