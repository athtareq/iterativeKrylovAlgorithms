function [P]=ProdDiam(Am,a,p)
  [m,l]=size(a);
  [n,mm]=size(Am);  
  P=a(1)*Am(:,1:p); 
  for i=2:l
    P=P+a(i)*Am(:,(i-1)*p+1:i*p);
  end
  
  
    
  
  