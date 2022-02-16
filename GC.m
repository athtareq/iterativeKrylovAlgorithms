function [x]=GC(A,b,m,x0)
  r0=b-A*x0;p0=r0;r=r0;
  eps=10^-6;
  while norm(r,2)>eps
    alpha=(r0'*r0)/((A*p0)'*p0);
    x=x0+alpha*p0;
    r=r0-alpha*(A*p0);
    beta=(r'*r)/(r0'*r0);
    p=r+beta*p0;
    r0=r;p0=p;
  end
  