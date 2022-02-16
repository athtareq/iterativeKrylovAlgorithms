function [P]=ProdDiam2(Am,M,p)
  [m,l]=size(M);
  [n,mm]=size(Am);  
  for i=1:l
    P(:,(i-1)*p+1:i*p)=ProdDiam(Am,M(:,i),p);
  end
  
