function [V,H]=nArnoldi_HH(A,m)
  [n,l]=size(A);
  r0=rand(n,1);
  z=r0;
  v=z+sign(z(1))*norm(z);
  v=v/norm(v);
  I=eye(n);
  P1=I-2*v*v';
  V(:,1)=P1*I(:,1);
  H=zeros(n,m);
  for j=1:m+1
    z=P1'*A*V(:,j);
    for i=1:j
      z(i)=0;
    end
    Pj=mat_Householder(z,j);
    Pj=P1*Pj;
    V(:,j+1)=P1*I(:,j+1);
    H(:,j)=Pj*z;
  end