function [v,h]=Arnoldi_HH(A,m)
  [n,l]=size(A);
  r0=rand(n,1);
  I=eye(n);
  Pjp=I;
  Pjm=I;
  z=r0;
  for j=1:m+1
    w=z/norm(z);
    if j>1 
    Pj=mat_Householder(w,j);
    Pjm=Pjm*Pj; %produit Pj dans ordre decroissant 
    Pjp=Pj*Pjp; %produit Pj dans ordre croissant   
    %pour simplifier les calcul prochains
    h(:,j-1)=Pj*z;
    end
    v(:,j)=Pjp*I(:,j);
    if j<=m
      z=Pjm*(A*v(:,j));
    end
  end
  
  
  