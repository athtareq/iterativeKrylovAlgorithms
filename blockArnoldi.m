function[Hm,Vm,x]=blockArnoldi(A,X0,m,B,p)
R0=B-A*X0;
[n,p]=size(R0);
Vm=zeros(n,p*(m+1));
Hm=zeros(p*(m+1),p*m);
I=eye(p);
solution=zeros(n,p);
V0=R0; Ipm=eye(p*m); II=eye(p*(m+1));
[Q,R]=qr(V0,0);
G=R;
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
Hmm=Hm(1:m*p,:); 
  [L,U]=lu(Hmm);
for i=1:p
    g=Ipm(:,1:p)*G*I(:,i);
    Z=L\g;
    Y=U\Z;
    X=Um*Y;
    solution(:,i)=X;
end

%sol=Hmm\(Ipm(:,1:p)*G);
x=solution;
erreur=norm(B-A*solution)
%Um*sol
end
