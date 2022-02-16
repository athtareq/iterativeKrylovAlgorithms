function [H,V]=ArnoldiGlobal(A,X0,m,B)
[n,n]=size(A);
R0=B-A*X0;
[n,p]=size(X0);
Vm=zeros(n,p*(m+1));
H=zeros((m+1),m);
V0=R0;Ipm=eye(p*m); II=eye(p*(m+1));
[Q,R]=qr(V0);
V=Q(:,1:p);V=R0/norm(R0);Vm(:,1:p)=V;
somme=zeros(n,p);
    for j=1:m
        W=A*Vm(:,(j-1)*p+1:j*p);
        for i=1:j
            H(i,j)=trace(Vm(:,(i-1)*p+1:i*p)'*W);
            w=W-H(i,j)*Vm(:,(i-1)*p+1:i*p);
        end
        H(j+1,j)=sqrt(trace(w'*w));
        if H(j+1,j)==0 exit(-1)
        end
        Vm(:,j*p+1:(j+1)*p)=w/H(j+1,j);        
    end
    Hm=H(1:m,:);
    Um=Vm(:,1:m);
    P=ProdDiam2(A,Um,p);
    err1=norm(TranspDiam(Um,P,m,p)-Hm)

end
    