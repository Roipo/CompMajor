function V = ComputeTutteParameterization(obj)
    Bi = obj.findOrientedBoundariesMatlab;
    Bi = Bi{1};
    t = linspace(0,2*pi,length(Bi)+1);
    Bv = [cos(t') sin(t')];
    Bv(end,:)=[];


    %% solve Tutte embedding for initial guess (assume nice angles of source mesh)
    L=obj.computeV2V;
    L=diag(sparse(sum(L,2)))-L;
    rhs = zeros(obj.Nv,2);
    rhs(Bi,:)=Bv;
    L(Bi,:)=0;
    L(Bi,Bi)=eye(numel(Bi));
    V = (L\rhs)';
end