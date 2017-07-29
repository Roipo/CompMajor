function [D1,D2] = computeSurfaceGradientMatrix(G,F1,F2)
% Computes the surface gradient matrix based on the frames F1,F2
%D1,D2 operators for surface gradients of size FXV - coordinates in triangle frame system
[Dx,Dy,Dz]=G.computeGradientMatrix;
if nargin<3
    [F1,F2] = G.getFramePerFace;
end
%F1=-F1;
% F1(1)Dx+F1(2)Dy+F1(3)Dz
%projecting gradient on triangle frame
% D1=diag(F1(1,:))*Dx+diag(F1(2,:))*Dy+diag(F1(3,:))*Dz;
% D2=diag(F2(1,:))*Dx+diag(F2(2,:))*Dy+diag(F2(3,:))*Dz;

D1=bsxfun(@times,F1(1,:)',Dx)+bsxfun(@times,F1(2,:)',Dy)+bsxfun(@times,F1(3,:)',Dz);
D2=bsxfun(@times,F2(1,:)',Dx)+bsxfun(@times,F2(2,:)',Dy)+bsxfun(@times,F2(3,:)',Dz);