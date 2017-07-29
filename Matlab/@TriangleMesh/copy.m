function copy(obj,m)
% COPY copies mesh object during constructor phase
%   input: obj - the TriangleMesh object
%          m   - the TriangleMesh object to copy
%
% Created by Nave Zelzer on may 22 2014.
obj.V 		 =  	m.V;
obj.Nv 		 = 		m.Nv;

obj.F 		 =  	m.F;
obj.Nf		 = 		m.Nf;

obj.E 		 =  	m.E;
obj.Ne 		 = 		m.Ne; 

obj.Material = 		m.Material;

end

