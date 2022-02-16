n=1000;
m=20;
A=4*eye(n)-diag(ones(n-1,1),-1)-diag(ones(n-1,1),+1);
 b1=A*ones(n,1);
 b2=2*(A*ones(n,1));
 b3=3*(A*ones(n,1));
 p=3;
 B=[b1,b2,b3];
 X0=zeros(n,p);
 R0=B-A*X0;
 %R0=B-A*X0;
 solution_exacte= [ones(n,1),2*ones(n,1),3*ones(n,1)]
[Hm,Vm,x]=blockArnoldi(A,X0,m,B,p);
for i=1:50
  [Hm,Vm,x]=blockArnoldi(A,X0,i,B,p);
  err(i)=norm(solution_exacte-x)
end
plot(err)