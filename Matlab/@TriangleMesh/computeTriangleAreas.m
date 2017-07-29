function TriArea = computeTriangleAreas( obj, EL )
% COMPUTESURFACEAREA computes the areas of all the triangles in the mesh
%   input: TriangleMesh object
%		   EL - edge lengths (optional)
%   outut:
%         TriArea - column vector of the areas of all the triangles in the
%                   mesh 
%
% Created by Nave Zelzer from code by Roi Poranne on march 17 2014.

% if EL (edge lengths) aren't given as input calculate them 
if exist('EL','var')==0
    EL=obj.computeTriangleEdgeLengths;
end

% calculating the areas using Heron's formula:
%      A2
%      /\
%   a /  \ b     S = 0.5*(a+b+c)
%    /____\      area = sqrt(S*(S-a)+S*(S-b)+S*(S-c))
%   A1 c  A3
L1=EL(:,1);
L2=EL(:,2);
L3=EL(:,3);
S=(L1+L2+L3)/2;
TriArea=sqrt(S.*(S-L1).*(S-L2).*(S-L3));
end

