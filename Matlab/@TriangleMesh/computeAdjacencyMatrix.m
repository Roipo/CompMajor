function A = computeAdjacencyMatrix(obj,type)
% COMPUTEADJACENCYMATRIX computes Adjacency Matrix for the mesh, thus the
%                         positive elements of each row i correspond to 
%                         the vertices that are incident to vi
%   input: the TriangleMesh object
%   outut: A the Adjacency Matrix
%
% Created by Roi Poranne.

E = obj.E;
Nv = obj.Nv;
I = [E(1,: ),E(2,:)];
J = [E(2,: ),E(1,:)];
if nargin==1
    type='combinatorial';
end

switch (type)
    case 'combinatorial'
        EL = ones(1,size(E,2));
    case 'distance'
        EL = computeEdgeLengths(obj);
end
A = sparse(I,J,[EL,EL],Nv,Nv);

end