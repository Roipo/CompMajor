function h = ComputeConeHessian(A1,A2,x,w)

Ax = [A1*x,A2*x];
f = sqrt(sum(Ax.^2,2));
wf=sparse(w./f);
wAtA=A1'*diag(wf)*A1+A2'*diag(wf)*A2;

AtAx=A1'*diag(sparse(A1*x))+A2'*diag(sparse(A2*x));
wf3=sparse(w./(f.^3));
wAtAx=AtAx*diag(wf3)*AtAx';

h=wAtA-wAtAx;