function h = ComputeConeHessianLoop(A1,A2,x)
% compute the hessians of the cones C_i = ||[A1(i,:);A2(i,:)]*x||
% H is a cell array where H{i} is the hessian of C_i
% This is the much slower loop version of ComputeConeHessianLoop

Ax = [A1*x,A2*x];
f = sqrt(sum(Ax.^2,2));
h=cell(1,size(A1,1));
for i=1:size(A1,1)
    A=[A1(i,:);A2(i,:)];
    AtA=A'*A;
    AtAx=AtA*sparse(x);
    h{i} = AtA/f(i)-(1/(f(i).^3))*(AtAx)*(AtAx');
end

