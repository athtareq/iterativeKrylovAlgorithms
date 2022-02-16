n=5000;m=15;In=eye(n);Im=eye(m);
A=3*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
b=A*ones(n,1);
x0=zeros(n,1);
for i=1:1:50
  [x]=GMRES(A,b,i,x0);
  err(i)=norm(b-A*x);
end
figure, plot(err);
title("Minimisation de l'erreur du residu par GMRES");
