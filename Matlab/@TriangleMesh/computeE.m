function E = computeE(obj)
% COMPUTEE computes E matrix. The edges are given as (i,j) where i<j, 
%           and are sorted lexicographically.
%   input: TriangleMesh object
%   output: E as described above
%
% Created by Nave Zelzer on may 11 2014.
F = obj.F;
% we pair up every triangle indices and create sparse matrix out of it
I = [F(1,:),F(2,:),F(3,:)];
J = [F(2,:),F(3,:),F(1,:)];
E = [I ; J];
% we sort every row such that (i,j) where i<j thus, E(i,1) < E(i,2)
% E = sort(E)';
Eswitch=E(1,:)>E(2,:);
E(:,Eswitch)=E([2;1],Eswitch);
E=E';
% sorting the rows lexicographically
E = sortrows(E);
% removing duplicates
E = unique(E,'rows')';