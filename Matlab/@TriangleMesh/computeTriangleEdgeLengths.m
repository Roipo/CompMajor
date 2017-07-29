function EL = computeTriangleEdgeLengths(obj)
% COMPUTEEDGELENGTHS Computes the triangle edge lengths of the mesh
%   input:  TriangleMesh object
%   output: matrix of the triangle edge lengths.
%
% Created by Roi Poranne, Updated by Nave Zelzer on march 17 2014.
V=obj.V';
F=obj.F';
% calculating side lengths of the triangles
%      A2
%      /\       side length 1 = sqrt(sum((A2 - A3)^2,2))
%   a /  \ b    side length 2 = sqrt(sum((A1 - A3)^2,2)) 
%    /____\     side length 3 = sqrt(sum((A1 - A2)^2,2)) 
%   A1 c  A3
L1 = sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2));
L2 = sqrt(sum((V(F(:,1),:)-V(F(:,3),:)).^2,2));
L3 = sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2));
EL = [L1, L2, L3];
end

