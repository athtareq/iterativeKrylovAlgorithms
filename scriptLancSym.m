n=1000;m=20;
v=rand(n,1);D=diag(v);u=rand(n,1);
u=u/norm(u);P=eye(n)-2*u*u';A=P*D*P';
Im=eye(m);
b=A*ones(n,1);
x0=zeros(n,1);
for i=1:50
  [T,V,x]=LanczosSym(A,b,i,x0);
  err(i)=norm(b-A*x);
end
figure,plot(err)

