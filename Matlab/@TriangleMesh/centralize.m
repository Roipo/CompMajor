function obj = centralize(obj,scale)
% CENTRALIZE Centres the mesh so its barycentre is at the origin
%   input: scale - a string that lets the user choose how to centre the
%                  mesh. there are three options:
%                       0) scale is not defined then just centre the mesh
%                       1) scale = 'ScaleArea' - scales the mesh so it will
%                                   have unit area.
%                       2) scale = 'BoundingBox' fits the mesh to [-1,1]^3
%                       3) scale = 'BoundingBox01' fits the mesh to [0,1]^3
%   outut:
%         pointer to the centred mesh
%
% Created by Roi Poranne, Updated by Nave Zelzer on march 3 2014.

if nargin==1
    scale=0;
end
Nv = obj.Nv;
% obj.V is a 3 by Nv matrix thus cxyz is a vector [xc, yc, zc] where xc
% is the barycentre of the x coordinates, yc for the y coordinates and zc
% for the z coordinates. the barycentre = (min+max)/2.
cxyz = (max(obj.V,[],2) + min(obj.V,[],2))/2;

% we subtract the barycentre from each vertex to centre it around the
% origin.
%     xc
% if  yc  is the barycentre we replicate it to Nv times and subtract 
%     zc
% it from V:
%
%              / v1x v2x ... vnx \     / xc  xc ... xc \
% centred V = |  v1y v2y ... vny  | - | yc  yc ...  yc | 
%              \ v1z v2z ... vnz /     \ zc  zc ... zc /
obj.V = obj.V - repmat(cxyz,1,Nv);

% we multiply the mesh by the square root of 1/A where A is the total area
% to transform to unit total area. 
if strcmp(scale,'ScaleArea')
    area = obj.computeSurfaceArea;
    obj.V = obj.V*sqrt(1/area);
end

% we normalize the vertices by the absolute maximum coordinate values to
% transform it to [-1,1]^3 instead of [-a,b]^3
if strcmp(scale,'BoundingBox')
    obj.V=obj.V./(max( abs(obj.V(:)) ));
end

% we normalize the vertices by the absolute maximum coordinate values 
% times 2 to transform it to [0,1]^3 instead of [-a,b]^3
if strcmp(scale,'BoundingBox01')
    obj.V=obj.V./(2*max(abs(obj.V(:))));
    obj.V=obj.V+0.5;
end