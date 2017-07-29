function constructFromObj(obj,objName)
% CONSTRUCTFROMOBJ copies mesh properties from obj file during constructor 
%                   phase
%   input: obj - the Mesh object
%          objName   - the name of the obj file
%
% Created by Nave Zelzer on may 22 2014.
s = TriangleMesh.objReaderMex(objName);
if s.isTriangleMesh == 0
    error('TriangleMesh object can take only triangle meshes');
end
obj.V = s.V';
obj.Nv = size(s.V,1);
obj.F = s.F';
obj.Nf = size(s.F,1);
obj.E = obj.computeE();
obj.Ne = size(obj.E,2); 

obj.Material.Vt = s.Vt';
obj.Material.Vn = s.Vn';
obj.Material.Ft = s.Ft';
obj.Material.Fn = s.Fn';

end

