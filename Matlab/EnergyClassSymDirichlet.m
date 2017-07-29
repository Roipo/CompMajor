classdef EnergyClassSymDirichlet < handle
    %Symmetric Dirichlet Energy gradient and quadratic proxies computations
    %class
    properties
        eps=1e-6 % eps for pushing PSD proxy to PD
    end
    properties(SetAccess = private)        
        Ai %elements areas
        AA %diagonal matrix of elements areas repeated 4 times

        D1 %derivative operator in frame direction 1 
        D2 %derivative operator in frame direction 2
        DD %nabla operator 

        x %u,v column stack
        
        fi %energy of terms of SOS objective        
        J
        dS
        ds

        %singular values of J (of the mapping)
        S
        s
        u1
        un
        v1
        vn
        
        %dirichlet proxy (constant w.r.t x)
        Hdirichlet
    end
    
    methods
       function obj=EnergyClassSymDirichlet()
           
       end
       
       function obj=Init(obj,D1,D2,Ai,V,F)
           %D1,D2 derivative operators in frame directions
           %Ai elements areas
           %vertices and faces: V is Nvx3, F is Nfx3
            obj.D1=D1;
            obj.D2=D2;
            obj.Ai=Ai;
            obj.DD=kron(eye(2),[obj.D1;obj.D2]);
            obj.AA=diag(sparse(repmat(obj.Ai,4,1)));
            obj.Hdirichlet = obj.DD'*obj.AA*obj.DD;           
            NewtonMex('Init',V,F);%used in Energy Class
        end
        
        function [f,g] = ComputeEnergy(obj)
            f=obj.Ai'*obj.ComputeEnergyPerElement;
            if  nargout>1
                g = (obj.Ai'*obj.ComputeGradientPerElement)';
            end
        end
        
        %method used for calculated scale for min sym-dirichlet energy
        function [fDirichlet,finDirichlet]=ComputeSymmetricDirichletParts(obj)
            fTerms=obj.ComputeEnergyTerms;
            sqrf = fTerms.^2;
            fDirichlet =    0.5*obj.Ai'*(sqrf(:,1) + sqrf(:,3));
            finDirichlet = 0.5*obj.Ai'*(sqrf(:,2) + sqrf(:,4));
        end      
        
        function UpdatePose(obj,x)
            %x is the current solution in the format [u,v] or [u ; v] (2 columns of column stack of parameterization)
            obj.x=x(:);
            obj.SetSingularValues;
        end
        
        %% Proxies
        function H=ProxySLIM(obj)
            Umat=[obj.u1,obj.un];
            Umat = transposeVec2x2(Umat);
            W=zeros(size(Umat,1),4);
            W(:,4) = sqrt(0.5*(obj.s-obj.s.^-3)./(obj.s-1));
            W(:,1) = sqrt(0.5*(obj.S-obj.S.^-3)./(obj.S-1));
            W(abs(obj.S-1)<1e-8,1)=1;
            W(abs(obj.s-1)<1e-8,4)=1;
            WW=sparse(multiplyVec2x2(multiplyVec2x2(Umat,W),transposeVec2x2(Umat)));
            
            WWWW=[diag([WW(:,1);WW(:,1)]),diag([WW(:,2);WW(:,2)]);diag([WW(:,3);WW(:,3)]),diag([WW(:,4);WW(:,4)])];
            
            WWWDD=WWWW*obj.DD;
            H=WWWDD'*obj.AA*WWWDD;
        end
        
        function H=ProxyCompMajor(obj)
            Js = [obj.dS;obj.ds];
            Hs = [(1+3*obj.S.^-4).*obj.Ai; (1+3*obj.s.^-4).*obj.Ai];
            Hggn = Js'*diag(sparse(Hs))*Js; %generalized gauss newton
            gS = obj.Ai.*(obj.S-obj.S.^-3);
            gs = obj.Ai.*(obj.s-obj.s.^-3);
            walpha = gS+gs;
            wbeta = gS-gs;
            walpha(walpha<0) = 0;

            a1=0.5*[obj.D1, obj.D2]; a2=0.5*[-obj.D2,obj.D1]; %similarity cone coefficients 
            b1=0.5*[obj.D1,-obj.D2]; b2=0.5*[ obj.D2,obj.D1]; %antisimilarity cone coefficients
            ha = ComputeConeHessian(a1,a2,obj.x,walpha);
            hb = ComputeConeHessian(b1,b2,obj.x,wbeta);
            
            H=Hggn+ha+hb;                       
        end
        
        % This is a reference implementation which is less vectorized but
        % might be easier to read
        function H=ProxyCompMajorLoop(obj)
            
            D1=obj.D1/2;
            D2=obj.D2/2;
            a1=[D1, D2]; a2=[-D2,D1];
            b1=[D1,-D2]; b2=[ D2,D1];
            ha = ComputeConeHessianLoop(a1,a2,obj.x);
            hb = ComputeConeHessianLoop(b1,b2,obj.x);
            
            s=obj.Energy.s; S=obj.Energy.S;
            [dS,ds]=obj.Energy.ComputeSingularDerivatives;
            Js = [dS;ds];
            Hs = [(1+3*S.^-4).*obj.Ai; (1+3*s.^-4).*obj.Ai];
            Hggn = Js'*diag(Hs)*Js;%generalized gauss newton
            gS = obj.Ai.*(S-S.^-3);
            gs = obj.Ai.*(s-s.^-3);

            walpha = gS+gs;
            wbeta  = gS-gs;
            walpha(walpha<0) = 0;

            Hab=zeros(size(Hggn));
            for i=1:length(ha)
                Hab=Hab+walpha(i)*ha{i};
                Hab=Hab+wbeta(i)*hb{i};
            end
            H=Hggn+Hab;
        end
        
        function H=ProxyNewton(obj)
            % bNoProjection - flag mitigating performing per face PSD projection
            %             bNoProjection = 0; % means not simple Newton and perform face projection.
            %             default is 0 (perform face projection)
            bNoProjection = 1;
            H = NewtonMex('Compute',obj.x,bNoProjection);            
        end
        
        function H=ProxyProjectedNewton(obj)
            H = NewtonMex('Compute',obj.x);
            %             For Hessians for each face, do this:
            %             [H,HH] = NewtonMex('Compute',x);
        end        
        
        function H=ProxyProjectedNewtonFull(obj)
           H = obj.ProxyNewton;
           [V,D] = eig(full(0.5*(H+H')));
           e = diag(D);
           e(e < obj.eps)= obj.eps;
           H = V*diag(e)*V';
        end        
    end
    
    methods(Access = private)
        function fe=ComputeEnergyPerElement(obj)
            obj.fi=obj.ComputeEnergyTerms;
            fe=0.5*sum(obj.fi.^2,2);
        end
        
        function f=ComputeEnergyTerms(obj)
            f=[obj.S,1./obj.S,obj.s,1./obj.s];
        end
        
        function ge=ComputeGradientPerElement(obj)
            obj.ComputeSingularDerivatives;
            Nf=numel(obj.S);
            
            dSE=[ones(Nf,1),-1./(obj.S.^2)]; dsE=[ones(Nf,1),-1./(obj.s.^2)];
            
            gS=(dSE(:,1).*obj.fi(:,1)+dSE(:,2).*obj.fi(:,2));
            gs=(dsE(:,1).*obj.fi(:,3)+dsE(:,2).*obj.fi(:,4));
            
            ge = diag(sparse(gS))*obj.dS+...
                diag(sparse(gs))*obj.ds;
        end
 
        function SetSingularValues(obj)
            Nv = length(obj.x)/2;
            U = obj.x(1:Nv);
            V = obj.x(Nv+1:end);
            
            % enetries of Jacobian of all elements (triangles) 
            a=obj.D1*U;
            b=obj.D2*U;
            c=obj.D1*V;
            d=obj.D2*V;
            
            %save for checking flips of elements
            obj.J=a.*d-b.*c;
            
            D1P=obj.D1*[U,V]/2; D2P=obj.D2*[U,V]/2;
            E=D1P(:,1)+D2P(:,2); F=D1P(:,1)-D2P(:,2);
            G=D2P(:,1)+D1P(:,2); H=-D2P(:,1)+D1P(:,2);
            Q=sqrt(E.^2+H.^2); R=sqrt(F.^2+G.^2);
            
            a1=atan2(G,F); a2=atan2(H,E);           
            theta=(a2-a1)/2; phi=(a2+a1)/2;
            
            %save singular values and vectors
            obj.s=Q-R;
            obj.S=Q+R;
            obj.u1=[  cos(phi)  sin(phi)];
            obj.un=[ -sin(phi)  cos(phi)];
            obj.v1=[ cos(theta) -sin(theta)];
            obj.vn=[ sin(theta)  cos(theta)];  
        end
        
        function ComputeSingularDerivatives(obj)
            b=diag(sparse(obj.v1(:,1)))*obj.D1+diag(sparse(obj.v1(:,2)))*obj.D2;
            c=diag(sparse(obj.vn(:,1)))*obj.D1+diag(sparse(obj.vn(:,2)))*obj.D2;
            obj.dS=[diag(sparse(obj.u1(:,1)))*b,diag(sparse(obj.u1(:,2)))*b];
            obj.ds=[diag(sparse(obj.un(:,1)))*c,diag(sparse(obj.un(:,2)))*c];
        end
    end
end

