n=1000;m=15;In=eye(n);Im=eye(m);
A=3*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
b=A*ones(n,1);
x0=zeros(n,1);
for i=1:50
  x=FOM(A,b,i,x0);
  err(i)=norm(ones(n,1)-x);
  y=GMRES(A,b,i,x0);
  err2(i)=norm(b-A*y,2);
  [T,V,z]=LanczosSym(A,b,i,x0);
  err3(i)=norm(b-A*z);  
end
figure, subplot(131), plot(err), title("FOM");
subplot(132), plot(err2),title("GMRES");
subplot(133), plot(err3),title("Lanczos");
B=2*In-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
b=B*ones(n,1);
for i=1:50
  x=FOM(B,b,i,x0);
  err(i)=norm(ones(n,1)-x);
  y=GMRES(B,b,i,x0);
  err2(i)=norm(b-A*y,2);
  [T,V,z]=LanczosSym(B,b,i,x0);
  err3(i)=norm(b-A*z);  
end
figure, subplot(131), plot(err), title("FOM");
subplot(132), plot(err2),title("GMRES");
subplot(133), plot(err3),title("Lanczos");

