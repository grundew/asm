classdef MultiLayerModel < handle
   
    properties
        
        % Target solid properties and thickness
        solid;
        thickness;
        
        % Fluid in front and back of target
        fluid;
                
    end
    
    methods
        
        function this = MultiLayerModel(arg, thickness)
            
            switch arg
                
                case 'watersteelwater'
                    this.solid = materials.MaterialFactory.produce('stainless steel');
                    f(1) = materials.MaterialFactory.produce('water');
                    f(2) = f(1);
                    this.fluid = f;
                    this.thickness = thickness;
                
                case 'aralditesteelwater'
                    % Araldite (fluid)
                    s = materials.Fluid();
                    s.name = 'araldite';
                    s.density = 2000;
                    s.v = 1800;
                    f(1) = s;
                    
                    % Water
                    s = materials.Fluid();
                    s.name = 'water';
                    s.density = 1000;
                    s.v = 1450;
                    f(2) = s;
                    
                    this.fluid = f;
                    % Steel used in Cervanka paper
                    s = materials.LinearElastic();
                    s.name = 'steel';
                    s.density = 7850;
                    s.v = 5900;
                    s.vShear = 3150;
                    this.solid = s;
                    this.thickness = thickness;
            end
            
            
        end
        
    end    
    
end