n=1000;m=15;In=eye(n);Im=eye(m);
A=3*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
b=A*ones(n,1);
x0=zeros(n,1);
x=FOM(A,b,m,x0);
x=restartedFOM(A,b,i,x0);
disp("Minimisation de l'erreur par Restarted FOM");
err=norm(ones(n,1)-x,2)
