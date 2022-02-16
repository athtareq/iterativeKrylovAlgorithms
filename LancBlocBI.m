function[V,W,erreur]=NVLancBlocBI(A,B,C,X0,m,p,v1,w1)
[n,n]=size(A);
VV=zeros(n,(m+1)*p);
WW=zeros(n,(m+1)*p);
R0=B-A*X0;
Rt=C-A'*X0;
z=v1'*w1; Imp=eye(m*p);
[L,U]=decompLU(z);
D=diag(diag(U));
%N=(inv(D))*U;
sqrtAbsD=diag(diag(sqrt(abs(D))));
signD=diag(diag(sign(D)));
delta=sqrtAbsD*U;
beta=L*sqrtAbsD*signD;
V=v1*(inv(delta));
W=w1*(inv(beta'));
VV(:,1:p)=V;
WW(:,1:p)=W;
V0=zeros(n,p);
W0=zeros(n,p);
beta=zeros(p,p); delta=zeros(p,p);
for j=1:m
    alpha=W'*(A*V);
    vt=A*V-V*alpha-V0*beta;
    wt=A'*W-W*(alpha')-W0*(delta');
    z=(wt')*vt;
    [L,U]=decompLU(z);
    D=diag(diag(U));
    %N=(inv(D))*U;
    sqrtAbsD=diag(diag(sqrt(abs(D))));
    signD=diag(diag(sign(D)));
    delta=(sqrtAbsD)*U;
    beta=L*sqrtAbsD*signD;
    V0=V;
    W0=W;
    V=vt*(inv(delta));
    W=wt*(inv(beta))';
    VV(:,(j-1)*p+1+p:j*p+p)=V;
    WW(:,(j-1)*p+1+p:j*p+p)=W;
end
wm=WW(:,1:m*p);
vm=VV(:,1:m*p);
erreur=norm(Imp-wm'*vm,'fro')
   