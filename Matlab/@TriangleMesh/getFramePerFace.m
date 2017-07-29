function [F1,F2,N] = getFramePerFace(obj)
V=obj.V;
F=obj.F;
E12=V(:,F(2,:))-V(:,F(1,:));
E13=V(:,F(3,:))-V(:,F(1,:));
N=cross(E12,E13);
N=N./repmat(sqrt(sum(N.^2,1)),3,1);
F1=E12./repmat(sqrt(sum(E12.^2,1)),3,1);
F2=cross(N,F1);
end

