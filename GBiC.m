function [P,x]=GBiC(A,b,m,c,x0)
  r0=b-A*x0;rt0=c-A'*x0;r=r0;rt=rt0;
  p0=r0;pt0=rt0;
  eps=10^-6;
  while norm(r,2)>eps
    alpha=(r0'*rt0)/((A*p0)'*pt0);
    x=x0+alpha*p0;
    r=r0-alpha*(A*p0);
    rt=rt0-alpha*(A'*pt0);
    beta=(r'*rt)/(r0'*rt0);
    p=r+beta*p0;
    pt=rt+beta*pt0;
    r0=r;rt0=rt;p0=p;pt0=pt;
  end
  