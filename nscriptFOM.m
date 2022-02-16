n=1000;m=15;In=eye(n);Im=eye(m);
A=3*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
b=A*ones(n,1);
x0=zeros(n,1);
x=FOM(A,b,m,x0);
for i=1:1:50
  x=FOM(A,b,i,x0);
  err(i)=norm(ones(n,1)-x,2);
end
figure, plot(err);
title("Minimisation de l'erreur par FOM");
