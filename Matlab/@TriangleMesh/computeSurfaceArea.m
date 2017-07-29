function Area = computeSurfaceArea(obj)
% COMPUTESURFACEAREA computes the total area of the mesh surface.
%   input: TriangleMesh object
%   output:
%         Area    - total area of mesh surface
%
% Created by Roi Poranne, Updated by Nave Zelzer on march 13 2014.

TriArea = obj.computeTriangleAreas;
Area=sum(TriArea);% total area