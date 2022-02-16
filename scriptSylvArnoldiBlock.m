n=1000;m=20;p=2;In=eye(n);Ip=eye(p,p);
A=2*In-diag(ones(n-1,1),-1)+diag(ones(n-1,1),1);
B=Ip;C=A*ones(n,p)+ones(n,p)*B;
X0=zeros(n,p);
[Hm,Vm,Xm]=SylvArnlodiBlock(A,X0,m,B,p,C);
sol=Xm;
norm(C-A*Xm-Xm*B)
