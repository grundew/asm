classdef MultiLayerModel < handle
    % Container class for fluid - solid - fluid models.
    %
    % Example 1:
    % mdl = MultiLayerModel('watersteelwater', 12.3e-3);
    %
    % Example 2:
    % 
    % fluid = materials.Fluid();
    % fluid.v = 1450;
    % fluid.density = 1000;
    % solid = materials.MaterialFactory.produce('stainless steel');
    % thickness = 12.3e-3;
    % mdl = MultiLayerModel(fluid, solid, fluid, thickness);
    
    properties
        
        % Target solid properties and thickness
        solid;
        thickness;
        
        % Fluid in front and back of target
        fluid;
                
    end
    
    methods
        
        function this = MultiLayerModel(varargin)
            % Constructor for MultiLayerModel class.
            %
            % obj = MultiLayerModel(arg, thickness);
            % obj = MultiLayerModel(fluidFront, solid, fluidBack, solidThickness);
            
            if nargin == 2            
                arg = varargin{1};
                this.thickness = varargin{2};
                
                switch arg
                
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
                this.fluid(1) = varargin{1};
                this.solid = varargin{2};
                this.fluid(2) = varargin{3};
                this.thickness = varargin{4};
            else
                error('DNV:WrongNumerOfArguments', 'Number of arguments must be either 2 or 3');
            end
            
        end
        
    end    
    
end