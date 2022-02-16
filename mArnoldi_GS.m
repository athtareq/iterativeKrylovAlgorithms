function [v,h]=mArnoldi_GS(A,m,r0)
  v(:,1)=r0/norm(r0);
  for j=1:m
    w=A*v(:,j);
    for i=1:j
      h(i,j)=w'*v(:,i);
      w=w-h(i,j)*v(:,i);
    end
    h(j+1,j)=norm(w);
    if h(j+1,j)==0
      break;
    else v(:,j+1)=w/h(j+1,j);
  end
end
  
      
      