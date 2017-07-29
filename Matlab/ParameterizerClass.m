classdef ParameterizerClass < matlab.mixin.Copyable
    %Parameterization class for disc triangle meshes
    %Does proxy iterations + line search
    
    properties
        M                           %triangle mesh class
        Energy                      %Energy class
        ProxyTypesList
        ComputeProxy                %Proxy method - calculating H using the compatible Proxy method from Energy class
        c_linesearch = 0.9          % factor for schaefer stepsize for first dengenerated triangle
        EnergyTypesList = {'Sym-Dirichlet','Neo-Hookean'};
        %
        EnergyType
    end
    
    properties(SetAccess = private)
        Ac                  %linear constraints
        %convergence parameters
        num_conv_iters = 5; %num of iterations in arow to fullfil ftolerance or xtolerance criteria
        ftol = 1e-6;
        xtol = 1e-10;
        fcounter
        xcounter
        
        %bar angle
        bar_angle = 110;
        
        ls_alpha = 0.2; % line search sufficinet decrease
        ls_beta = 0.5;  % line search step factor
    end
    
    methods
        function obj=ParameterizerClass
            obj.SetEnergyType('Sym-Dirichlet')
        end
        
        function SetEnergyType(obj,Type)
            obj.EnergyType = Type;         
            if ~isempty(obj.Energy)
                obj.Energy.delete;
            end
            switch obj.EnergyType
                case 'Sym-Dirichlet'
                    obj.ProxyTypesList = {'Composite Majorization','SLIM','Projected Newton', 'Full projected Newton','Newton'};
                    obj.Energy = EnergyClassSymDirichlet;
                otherwise
                    error('incorrect energy type')
            end
            obj.SetProxyType(obj.ProxyTypesList{1})
        end
        
        function K = SetMesh(obj,M) % performs initialize too  
            obj.M = M;
            [D1,D2] = obj.M.computeSurfaceGradientMatrix;
            Ai=obj.M.computeTriangleAreas;
            obj.Energy.Init(D1,D2,Ai,obj.M.V',obj.M.F');
            K = obj.Initialize(D1,D2,Ai);
            obj.ResetOptimization;
        end
        
        function ResetOptimization(obj)
            obj.fcounter = 0;
            obj.xcounter = 0;
        end
        
        function setAc(obj)
            %% 3 coordinates
            %             s = [1,1,1];
            %             i = [1,2,3];
            %             j = [1,obj.M.Nv+1,obj.M.Nv+2];
            
            %% 1 point
            i = [1,2];
            j = [1,1+obj.M.Nv];
            s = [1,1];
            %% all boundary
            %             Bi = cell2mat(K.findOrientedBoundariesMatlab);
            %             s = ones(1,length(Bi)*2);
            %             i = 1:length(Bi)*2;
            %             j = [Bi,Bi+obj.M.Nv];
            
            %% sum zero
            %              s = ones(1,obj.M.Nv*2);
            %              i = [ones(1,obj.M.Nv),2*ones(1,obj.M.Nv)];
            %              j = [1:obj.M.Nv,(1:obj.M.Nv)+obj.M.Nv];
            obj.Ac = sparse(i,j,s,max(i),2*obj.M.Nv);
        end
        
        function SetProxyType(obj,type)
            e = obj.Energy;
            switch type
                case 'Composite Majorization'
                    obj.ComputeProxy = @e.ProxyCompMajor;
                case 'SLIM'
                    obj.ComputeProxy = @e.ProxySLIM;
                case 'Projected Newton'
                    obj.ComputeProxy = @e.ProxyProjectedNewton;
                case 'Full projected Newton'
                    obj.ComputeProxy = @e.ProxyProjectedNewtonFull;
                case 'Newton'
                    obj.ComputeProxy = @e.ProxyNewton;
            end
        end
        
        function [K,f,g]=Initialize(obj,D1,D2,Ai)
            switch obj.EnergyType
                case 'Sym-Dirichlet'
                    %initialize with Tutte embedding into a circle
                    K=TriangleMesh('VF',obj.M.ComputeTutteParameterization,obj.M.F);
                    obj.setAc;
                    [f,g]=ComputeEnergy(obj,reshape(K.V',[],1));
                    %calculated scale for min sym-dirichlet energy
                    [fDirichlet,finvDirichlet]=obj.Energy.ComputeSymmetricDirichletParts();
                    scale = (finvDirichlet/fDirichlet)^ 0.25;
                    K.V = scale*K.V;
                otherwise
                    error('incorrect energy type')
            end
        end
        
        function [x,f,g,dx,stepSize]=DoIteration(obj,x)
            %x is the current solution in the format [u,v] or [u ; v] (2 columns of column stack of parameterization)
            %output x should be of 2 columns [u,v]
            if size(x,2)==2
                x=x(:);
                reshapeFlag=true;
            else
                reshapeFlag=false;
            end
            %here x is column stack of parameterization
            [f,g,H]=ComputeEnergy(obj,x);
            Nconstraints =size(obj.Ac,1);
            KKT = [H, obj.Ac'; obj.Ac, zeros(Nconstraints) ];
            %             condest(KKT)
            rhs = [-g; zeros(Nconstraints,1)];
            sol = KKT\rhs;
            p = sol(1:2*obj.M.Nv);
            
            [t,stepSize,f] = obj.DoLineSearch(x,p,f,g);
            x=x+t*p;
            dx=norm(t*p);
            
            %make sure output x is 2 columns
            if reshapeFlag
                x=reshape(x,[],2);
            end
        end
        
        function [t,stepSize,f_new] = DoLineSearch(obj,x,p,f,g)
            %INPUT:
            %x - curr solution
            %p - step direction
            %f - curr objective value
            %g - gradient
            %OUTPUT
            %t - final step size
            %stepSize - step to first degenerate triangle (max value is 1), by bijective parameterization with free boundary by Smith & Schaefer
            %fnew - obj value after performing step
            
            stepSize = computeInjectiveStepSize(obj.M.F',x,p,1e-12);
            t = min(obj.c_linesearch*stepSize,1);
            
            
            alpha_g_p = obj.ls_alpha*g'*p;
            
            [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
            while linesearch_cond_lhs>linesearch_cond_rhs
                t = obj.ls_beta*t;
                [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
            end
            f_new = linesearch_cond_lhs;
            
            function [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond
                linesearch_cond_lhs = obj.ComputeEnergy(x + t*p);
                linesearch_cond_rhs = f + t*alpha_g_p;
            end
            
        end
        
        function [f,g,H]=ComputeEnergy(obj,x)
            %INPUT:
            %x - curr solution
            %OUTPUT:
            %f - curr objective value
            %g - gradient
            %H - quadratic Proxy Hessian to be used in KKT
            
            obj.Energy.UpdatePose(x);
            
            if  nargout==1
                f = obj.Energy.ComputeEnergy;
                return;
            end
            if nargout > 1
                [f,g]=obj.Energy.ComputeEnergy;
            end
            if nargout > 2
                H=obj.ComputeProxy();
            end
        end
        
        function status = OptimizationConverged(obj, fcur , fprev, xcur, xprev)
            status = 'Not Converged';
            if (abs(fcur - fprev) < obj.ftol*(fcur + 1))
                if (obj.fcounter >= obj.num_conv_iters)
                    status = 'Change in energy < tol';
                    return;
                else
                    obj.fcounter = obj.fcounter + 1;
                end
            else
                obj.fcounter = 0;
            end
            
            if (norm(xcur - xprev) < obj.xtol*(norm(xcur) + 1))
                if (obj.xcounter >= obj.num_conv_iters)
                    status = 'Change in X < tol';
                    return;
                else
                    obj.xcounter = obj.xcounter + 1;
                end
            else
                obj.xcounter = 0;
            end
        end
        
        function delete(obj)
            clear mex;
        end
    end
end

