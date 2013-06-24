classdef MultiLayerModel < handle
    % Container class for fluid - solid - fluid models.
    %
    % Example 1:
    % d_Tx = 10e-2;
    % d_Rx = 5e-2;
    % a_Tx = 6e-3;
    % a_Rx = 9e-3;
    % mld = MultiLayerModel('watersteelwater', 12.3e-3, d_Tx, d_Rx, a_Tx, a_Rx);
    %
    % Example 2:
    % 
    % fluid = materials.Fluid();
    % fluid.v = 1450;
    % fluid.density = 1000;
    % solid = materials.MaterialFactory.produce('stainless steel');
    % thickness = 12.3e-3;
    % d_Tx = 10e-2;
    % d_Rx = 8e-2;
    % a_Tx = 6e-3;
    % a_Rx = 9e-3;
    % alphaLambda = 9.3e-3;
    % mld = MultiLayerModel(fluid, solid, fluid, thickness,...
    %         d_Tx, d_Rx, a_Tx, a_Rx, alphaLambda);
    
    properties
        
        % Target solid properties and thickness
        solid;
        thickness;
        
        % Fluid in front and back of target (vector)
        fluid;

        % Alpha (damping factor)
        alphaLambda;

        % Distance from transmitter to plate
        d_Tx;
        % Distance from receiver to plate
        d_Rx;
        % Radius of transmitter
        a_Tx;
        % Radius of receiver
        a_Rx;
        
    end
    
    methods
        
        function this = MultiLayerModel(varargin)
            % Constructor for MultiLayerModel class.
            %
            % obj = MultiLayerModel(arg, thickness, d_Tx, d_Rx, a_Tx, a_Rx, alphaLambda)
            % obj = MultiLayerModel(fluidFront, solid, fluidBack, thickness, d_Tx, d_Rx)
            % obj = MultiLayerModel(fluidFront, solid, fluidBack, thickness);
            % obj = MultiLayerModel(fluidFront, solid, fluidBack, thickness, d_Tx, d_Rx, alphaLambda)
            
            if nargin == 7
                % obj = MultiLayerModel(arg, thickness, d_Tx, d_Rx)
                arg = varargin{1};
                this.thickness = varargin{2};
                this.d_Tx = varargin{3};
                this.d_Rx = varargin{4};
                switch arg
                
                    case 'naturalgassteel'
                        this.solid = materials.MaterialFactory.produce('stainless steel');
                        f(1) = materials.MaterialFactory.produce('natural gas');
                        f(2) = f(1);
                        this.fluid = f;

                    case 'naturalgassteelbitumen'
                        this.solid = materials.MaterialFactory.produce('stainless steel');
                        f(1) = materials.MaterialFactory.produce('natural gas');
                        f(2) = materials.MaterialFactory.produce('bitumen');
                        this.fluid = f;
                    
                    case 'watersteelwater'
                        this.solid = materials.MaterialFactory.produce('stainless steel');
                        f(1) = materials.MaterialFactory.produce('water');
                        f(2) = f(1);
                        this.fluid = f;
                        
                    case 'aralditesteelwater'
                        % Reference:
                        % Cervanka and Challande: A new efficient algorithm to
                        % compute the exact reflection and transmission factors
                        % for plane wase in layered absorbing media (liquid and
                        % solids)
                        
                        % Araldite (fluid)
                        s = materials.Fluid();
                        s.name = 'araldite';
                        s.density = 2000;
                        s.v = 1800;
                        s.alpha = 1000;
                        f(1) = s;
                        
                        % Water
                        s = materials.Fluid();
                        s.name = 'water';
                        s.density = 1000;
                        s.v = 1450;
                        s.alpha = 1;
                        f(2) = s;
                        
                        this.fluid = f;
                        % Steel used in Cervanka paper
                        s = materials.LinearElastic();
                        s.name = 'steel';
                        s.density = 7850;
                        s.v = 5900;
                        s.vShear = 3150;
                        s.alpha_L = 1000;
                        s.alpha_S = 1000;
                        this.solid = s;
                end
                
            elseif nargin == 4
                % obj = MultiLayerModel(fluidFront, solid, fluidBack, thickness);
                f(1) = varargin{1};
                f(2) = varargin{3};
                this.fluid = f;
                this.solid = varargin{2};
                this.thickness = varargin{4};
            elseif nargin == 9
                % obj = MultiLayerModel(fluidFront, solid, fluidBack,...
                %             solidThickness, d_Tx, d_Rx, a_Tx, a_Rx, alphaLambda)
                f(1) = varargin{1};
                f(2) = varargin{3};
                this.fluid = f;
                this.solid = varargin{2};
                this.thickness = varargin{4};
                this.d_Tx = varargin{5};
                this.d_Rx = varargin{6};
                this.a_Tx = varargin{7};
                this.a_Rx = varargin{8};
                this.alphaLambda = varargin{9};
            else
                error('HW:WrongNumerOfArguments',...
                    'Number of arguments must be either 6, 8 or 9');
            end
            
        end
        
    end
    
end