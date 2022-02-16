function [Hm,Vm,Xm]=SylvArnlodiBlock(A,X0,m,B,p,C)
R0=C-A*X0-X0*B;
[n,p]=size(R0);
Vm=zeros(n,p*(m+1));
Hm=zeros(p*(m+1),p*m);
I=eye(p);
solution=zeros(n,p);
V0=R0; Ipm=eye(p*m); II=eye(p*(m+1));
[Q,R]=qr(V0,0);
C1=R;
V=Q(:,1:p); Vm(:,1:p)=V;
for j=1:m
    w=A*V; EE=Ipm(:,(j-1)*p+1:j*p);
    for i=1:j
        E=II(:,(i-1)*p+1:i*p);
        V=Vm(:,(i-1)*p+1:i*p);
        H=V'*w;
        Hm=Hm+E*H*EE';
        w=w-V*H;
    end
    [Q,R]=qr(w,0);
    V=Q;
    H=R;
    E=II(:,j*p+1:(j+1)*p);
    Hm=Hm+E*R*EE';
    Vm(:,j*p+1:(j+1)*p)=V;
end
Um=Vm(:,1:m*p);
Hmm=Hm(1:m*p,1:m*p);
%erreur=norm(A*Um-Um*Hm-Vm(:,(m)*p+1:(m+1)*p)*Hm(m*p:(m+1)*p,(m-1)*p+1:m*p),'fro')
erreur1=norm(Um'*Um-Ipm)
erreur2=norm(A*Um-Um*Hm(1:m*p,:)-Vm(:,m*p+1:(m+1)*p)*Hm(m*p+1:(m+1)*p,(m-1)*p+1:m*p)*Ipm(:,(m-1)*p+1:m*p)')
%erreur=norm(Um'*Um-Ipm,'fro')
erreur3=norm(Um'*A*Um-Hmm)
Ym=sylvester(Hmm,B,Ipm(:,1:p)*C1);
Xm=X0+Um*Ym;
