classdef TriangleMesh < handle
    
    %MESH General purpose TRIANGLE mesh class for geometrical computations
    %   Multi-purpose mesh class for topological and geometrical mesh processing
    
    properties
        F  % 3 by Nf matrix. every column j contains the indices of the vertices on the j'th face
        Nf % number of faces
        V  % 3 by Nv matrix. every column i contains coordinates of the i'th vertex
        Nv % number of vertices
        E  % 2 by Ne matrix. every column j contains the indices of the vertices on the j'th edge
        Ne % number of edges
        Material % collects all the normal and texture data of the mesh
    end
    
    
    properties(Access = private)
        
    end
    
    % public methods
    methods
        function obj = TriangleMesh(varargin)
            obj.init();
            % copy constructor
            if (length(varargin)==1) && isa(varargin{1},'TriangleMesh')
                obj.copy(varargin{1});
            elseif length(varargin)>=2
                switch(varargin{1})
                    case 'obj'
                        s=obj.objReader(varargin{2});
                        if s.isTriangleMesh == 0
                            error('TriangleMesh object can take only triangle meshes');
                        end
                        obj.V = s.V';
                        obj.F = s.F';
                        
                        obj.Material.Vt = s.Vt';
                        obj.Material.Vn = s.Vn';
                        obj.Material.Ft = s.Ft';
                        obj.Material.Fn = s.Fn';
                    case 'off'
                        [V,F,Fs]=obj.offReader(varargin{2});
                        if any(Fs~=3)
                            error('TriangleMesh object can take only triangle meshes');
                        end
                        obj.V=V;
                        obj.F=F;
                    case 'VF'
                        obj.V=varargin{2};
                        obj.F=varargin{3};
                    
                
                        
                end
                obj.Nv = size(obj.V,2);
                obj.Nf = size(obj.F,2);
                obj.E = obj.computeE;
                obj.Ne = size(obj.E,2);

            end
        end
    end
    % public static methods implementations
    methods(Static)
        o = objReader(filename)
        [V,F,Fs] = offReader(filename)
    end % end of public static methods
    
    % private methods
    methods(Access = private)
        % mesh constructor initialization functions
        init(obj)
        copy(obj,m)
        constructFromObj(obj,objName)
        % mesh auxiliary topology and data structures functions
        E = computeE(obj)

    end % end of private methods
    
    % private static methods
    methods(Static, Access = private)
        % mex file wrappers
        ps = objReaderMex(objName)
    end % end of private static methods
    
end % end of class definition


