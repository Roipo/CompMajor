function Fn=computeFaceNormals(obj)
% COMPUTEFACENORMALS computes normals for each of the faces on the mesh
%   input: the TriangleMesh object
%   outut: 3 by nf matrix of normal vectors
%
% Created by Roi Poranne, Updated by Nave Zelzer on march 30 2014.

% old code for further testing
% F=obj.F;
% V=obj.V;
% nf=size(F,2);
% Fn=zeros(3,nf);
% if size(F,1)~=3
%     error('Not a triangular mesh!');
% end
% for i=1:nf
%     e1=V(:,F(2,i))-V(:,F(1,i));
%     e2=V(:,F(3,i))-V(:,F(1,i));
%     Fn(:,i)=cross(e1,e2);
%     Fn(:,i)=Fn(:,i)./sqrt(sum(Fn(:,i).^2));
% end

%       __         __
%  e1 = AB,   e2 = AC for each of the triangles in F
%       A
%       /\     
%      /  \
%     /____\
%    B      C
F=obj.F;
V=obj.V;
e1=V(:,F(2,:))-V(:,F(1,:));
e2=V(:,F(3,:))-V(:,F(1,:));
Fn = cross(e1,e2);
% we calculate the norm for each normal (row of norms) and replicate it
% three times so we could divide each normal coordinate by it.
Fn = Fn./repmat(sqrt(dot(Fn,Fn)),3,1);