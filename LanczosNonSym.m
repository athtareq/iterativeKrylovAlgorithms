function [x,xt]=LanczosNonSym(A,b,c,m,x0)
  [n,n]=size(A);Im=eye(m);
  r0=b-A*x0;  r0p=c-A'*x0;%r0p c'est r0' (prime)
  d=sqrt(abs(r0'*r0p));  B=(r0'*r0p)/d;
  v0=zeros(n,1);w0=zeros(n,1);
  v=r0/d;w=r0p/B;
  V=zeros(n,m+1);W=zeros(n,m+1);T=zeros(m+1,m);
  V(:,1)=v;W(:,1)=w;
  beta=0;delta=0;
  for k=1:m
    if k>1
      T(k-1,k)=beta;
      S(k-1,k)=delta;
    end
    alpha=w'*(A*v);
    t=A*v-beta*v0-alpha*v;
    z=A'*w-delta*w0-alpha*w;
    delta=sqrt(abs(z'*t));
    beta=(z'*t)/delta;
    T(k,k)=alpha;
    S(k,k)=alpha;
    v0=v;w0=w;
    v=t/delta;
    w=z/beta;
    V(:,k+1)=v;W(:,k+1)=w;
    T(k+1,k)=delta;S(k+1,k)=beta;
  end 
  %test des erreurs
  err1=norm(A*V(:,1:m)-V*T);
  err2=norm(A'*W(:,1:m)-W*S);
  sca=W(:,1:m)'*V(:,1:m);
  err_orthWV=norm(sca-eye(m));
  err_WtAV=norm(W(:,1:m)'*A*V(:,1:m)-T(1:m,:));
  %resolution
Tm=T(1:m,:);
Sm=S(1:m,:);
%trouver y et z
[L,U]=lu(Tm);
y=L\(B*Im(:,1));
ym=U\y;
x=x0+V(:,1:m)*ym;
[L,U]=lu(Sm);
z=L\(d*Im(:,1));
zm=U\z;
xt=x0+W(:,1:m)*zm;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

    