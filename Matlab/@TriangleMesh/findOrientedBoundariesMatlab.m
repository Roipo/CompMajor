function [OBV,OBE]  = findOrientedBoundariesMatlab( obj )
% FINDORIENTEDBOUNDARIESMATLAB returns a cell array of oriented boundary 
%  vertices (a cell for every boundary) and a cell array for boundary edges
%  a cell for every boundary.
%   input: TriangleMesh object
%   output: OBV cell array with a cell for every boundary vertices
%   output: OBE cell array with a cell for every boundary edges
%
% Created by Nave Zelzer on may 23 2014.
[~,BE] = obj.findBoundaries();
if isempty(BE)
    OBV = {};
    OBE = {};
    return;
end
OBV = {[]};
OBE = {zeros(2,0)};


bi = 1;
ind = 1;
while true
    OBV{bi}(end+1) = BE(1,ind);
    OBE{bi}(:,end+1) = BE(:,ind);
	ind1 = find(BE(1,:)==BE(2,ind));
    BE(1,ind)=0;
    if isempty(ind1) % close a boundary loop
        if ~any(BE(1,:))
            break;
        end
        OBV{end+1}=[];
        OBE{end+1} = zeros(2,0);
        bi = bi+1;
        ind = find(BE(1,:),1);
        continue;
    end
    ind=ind1;
end