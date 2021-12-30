%
%  Hx = H(x): computes the entropy of the integer valued matrix x
%

function Hx=H(x)
    x=reshape(x',1,[]);
    m=min(x);
    M=max(x);
    l=length(x);
    
    N=zeros(1,M-m+1);
    for i=1:l
        N(x(i)-m+1)=N(x(i)-m+1)+1;
    end
    
    nz=find(N);
    Hx=-N(nz)/l*(log2(N(nz)/l))';
end

