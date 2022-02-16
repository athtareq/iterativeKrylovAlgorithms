function[x]=restartedFOM(A,b,m,x0)
r0=b-A*x0;
beta=norm(r0);
Im=eye(m);
[n,n]=size(A);
tolerance = 10*e^-9;
k=1;
[V,H]=mArnoldi_GS(A,m,r0);
Hm=H(1:m,1:m);
Vm=V(:,1:m);
[L,U]=lu(Hm);
while beta>tolerance 
    %erreur=abs(y(m)*H(m+1,m));
    w=L\(beta*Im(:,1));
    y=U\w;
    x=x0+Vm*y;
    r0=r0-A*Vm*y;
    beta=norm(r0,2);
    x0=x;
    k=k+1;
end