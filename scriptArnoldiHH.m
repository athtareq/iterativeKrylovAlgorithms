n=5000;m=15;In=eye(n);Im=eye(m);
A=3*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
[v,h]=Arnoldi_HH(A,m);
Vm=v(:,1:m);
Hm=h(1:m,1:m);
Vm1=v(:,m+1)*Im(:,m)';
norm(A*Vm-Vm*Hm-h(m+1,m)*Vm1,'fro')