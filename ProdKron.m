function [C]=ProdKron(A,B)
  [n,p]=size(A);[m,q]=size(B);
  C=zeros(n*m,p*q);  
  for i=1:n
    for j=1:p
      C((i-1)*m+1:i*m,(j-1)*q+1:j*q)=A(i,j)*B;
    end
  end
  