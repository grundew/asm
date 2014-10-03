classdef TestWater < matlab.unittest.TestCase
    %   Test cases for startAsm

    properties
        
        fixture;
        testFigure;
        params;
        
        % Parameters to check
        nfft;
        
        % Results
        V
        
        % Computation time
        cputime
        
    end
    
    methods (TestClassSetup)
       
        function setupParams(this)
            this.nfft = 500;
            f = linspace(100e3, 300e3, this.nfft);
            p = parseAsmInput('cp', 1500,...
                'rho_fluid', 1000, 'f', f, 'savemat', false);
            this.params = p;
            tic
            this.V = computeAsmIntegral(@integrandFluidSolidFluid_planepiston, p);
            this.cputime = toc
        end
        
    end
    
    methods (Test)
        
        function testParameters(this)
            p = this.params;
            
            this.assertLength(this.V, this.nfft);
            this.assertSize(p.f, [1, this.nfft]);
        end
        
        function visualInspection(this)
            
            p = this.params;
            figure
            plot(p.f, abs(this.V))
            figure
            plot(p.f, unwrap(angle(this.V)))

        end
         
        function testValue(this)
            
            res = load('testValue_water.mat');
            this.verifyEqual(this.V, res.pt)
            
        end
        
        function testCpuTime(this)
            cpuupperlimit = 7;
            this.verifyLessThanOrEqual(this.cputime, cpuupperlimit);
        end
        
    end
    
end

