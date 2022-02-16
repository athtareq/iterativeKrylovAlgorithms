n=1000;m=15;In=eye(n);Im=eye(m);
A = 2*In-diag(ones(n - 1, 1), -1) + diag(ones(n-1, 1), +1);
b=A*ones(n,1);
c=A'*ones(n,1);
x0=zeros(n,1);
for i=10:100
  [x,xt]=LanczosNonSym(A,b,c,i,x0);
  err_sys(i)=norm(b-A*x);
  err_dual(i)=norm(c-A'*xt);  
end
figure, plot(err_sys)
figure, plot(err_dual) 





