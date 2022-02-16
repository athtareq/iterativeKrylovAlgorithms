function [M]=TranspDiam(W,V,m,p)
  [n,z]=size(W);
  M=zeros(m,m);
  for i=1:m
    for j=1:m
      M(i,j)=trace(W(:,(i-1)*p+1:i*p)'*V(:,(j-1)*p+1:j*p));
    end
  end
