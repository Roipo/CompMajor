function V2V = computeV2V(obj)
% COMPUTEADJACENCYMATRIX syntactic sugar for computeAdjacencyMatrix with
%                         more intuative name
%   input: the Mesh object
%   outut: V2V the Adjacency Matrix
%
% Created by Nave Zelzer on may 22 2014.
V2V = obj.computeAdjacencyMatrix();
