function [T,V,x]=LanczosSym(A,b,m,x0)
  [n,n]=size(A);
  T=zeros(m+1,m);V=zeros(n,m+1);
  r0=b-A*x0;B=norm(r0);v=r0/B;v0=zeros(n,1);
  beta=0;V(:,1)=v;
  Im=eye(m);  
  for k=1:m
    alpha=v'*A*v;
    T(k,k)=alpha;
    w=A*v-beta*v0-alpha*v;
    beta=norm(w);v0=v;
    v=w/beta;
    V(:,k+1)=v;T(k+1,k)=beta;
  end
beta=norm(b-A*x0);
[L,U]=lu(T(1:m,:));
y=L\(beta*Im(:,1));
ym=U\y;
x=x0+V(:,1:m)*ym;
  
    
