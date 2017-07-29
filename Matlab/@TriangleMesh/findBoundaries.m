function [BV,BE] = findBoundaries(obj)
% find the boundary vertices of G
% out: BV: indices of the boundary vertices
% out: BE: indices of the boundary edges
% Created by Nave Zelzer on may 22 2014.
F=obj.F;
Ic = [F(1,:),F(2,:),F(3,:)];
J = [F(2,:),F(3,:),F(1,:)];
E = [Ic ; J];
% we sort every row such that (i,j) where i<j thus, E(i,1) < E(i,2)
Es = sort(E);

% returns edge index for every edge in EFA, thus if we have duplicates the
% indices will have duplicates too.
[~, Ia, Ic] = unique(Es(1:2,:)', 'rows');
% count how many of each edge
binc = histcounts(Ic,(1:max(Ic)+1)-0.5);
% if only one edge in the count we have a boundary.
BE = E(:,Ia(binc==1));
% get only the boundary vertices
BV = unique(BE(:))';