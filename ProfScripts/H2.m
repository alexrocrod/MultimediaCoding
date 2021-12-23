%
%  H2x = H2(x): computes the per symbol second order entropy H(X1,X2)/2 
%              of the integer valued matrix x (elemnts are taken row-wise)
%

function H2x=H2(x)
    x=reshape(x',1,[]);
    m=min(x);
    M=max(x);
    l=length(x);
    
    N=zeros(M-m+1,M-m+1);
    for i=1:l-1
      N(x(i)-m+1,x(i+1)-m+1)=N(x(i)-m+1,x(i+1)-m+1)+1;
    end
    
    nz=find(N);
    H2x=-(1/2)*sum((N(nz)/(l-1)).*(log2(N(nz)/(l-1))));
end
